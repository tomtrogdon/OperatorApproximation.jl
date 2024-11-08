#LagSeries calls Lag to collect Laguerre polynomials
function Lag(n::Int64,x::Float64) # evaluate Laguerre polynomials at x, α = 1
    c = zeros(Float64,n+2)
    c[2] = exp.(-x/2)
    k = 0
    for i = 3:n+2 #k = 1 is the zeroth order, gives (k-1)th order
      c[i] = (2*k+2-x)*c[i-1]/(k+1) - c[i-2]
      k = k + 1
    end
    return c
end

#Res we call calls this LagSeries; effectively returns r_{j,α}(z) from AKNS paper
function LagSeries(j::Int64,z::Float64,x::AbstractVector{T}) where T # need to investigate stability
    c = Lag(j,z) # c[2] gives L_0^(1), c[i] gives L_{i-2}^(1)
    out = zeros(Complex{Float64},length(x),j)
    out[1:end,1] = x*c[2]
    @inbounds for i = 2:j
      out[1:end,i] = (1.0 .+ x).*out[1:end,i-1] .+ x*(c[i+1]-c[i])
    end
  
    return out
end

function LagSeries(z::Complex{Float64},x::Vector{Complex{Float64}},cfs::Vector{Complex{Float64}}) # need to investigate stability
    j = length(cfs)
    c0 = 0.
    c1 = exp(-z/2) # c[2]
    out = x*c1
    k = 0
    # pm = -1
    # sum = (pm*cfs[1])*out
    sum = cfs[1]*out
    @inbounds for i = 2:j
    #   pm = pm*(-1)
      c2 = (2*k+2-z)*c1/(k+1) - c0 # c2 = c[i + 1]
      out .*= (1.0 .+ x)
      out .+= x*(c2-c1)
    #   axpy!(pm*cfs[i],out,sum)
    axpy!(cfs[i],out,sum) #BLAS version of sum += cfs[i]*out
      k += 1
      c0 = c1
      c1 = c2
    end
    return sum
end

function Res(j::Int64,α::Float64,z::Vector{Complex{Float64}})
    x = (-2im*sign(j))./(z.+1im*sign(j))
    y = 2*abs(α)
    return -LagSeries(abs(j),y,x)
end

#CauchyPNO calls this one
function Res(j::Int64,α::Int64,z::GridValues)
    x = (-2im*sign(j))./(z.+1im*sign(j))
    y = 2*abs(α)
    return -LagSeries(abs(j),y,x)
end
  
function Res(j::Integer,α::Complex{Float64},z::Vector{Float64},cfs::Vector{Complex{Float64}})
    x = (-2im*sign(j))./(z.+1im*sign(j))
    y = 2*abs(α)
    return -LagSeries(y,x,cfs)
end

#effectively returns -M_{+1}(k) in α<0 case and M_{-1}(k) in α>0 case from AKNS paper
function CauchyPNO(n,m,α,z)
    α = convert(Float64,α)
    if α < 0.
        return -Res(m,α,z(n))
    else
        return Res(-m,α,z(n))
    end
end

function CauchyConstantMatP(i,j)
    if i == j
        if j <= 0
            return 0
        else 
            return 1
        end
    elseif (i == 0) & (j > 0)
        return -1
    else
        return 0
    end
end

function CauchyConstantMatM(i,j)
    (i == j ? -1 : 0) + CauchyConstantMatP(i,j) #C- = C+ - I
end

function BuildOperatorBlock(n,m,α,gridPts)
    A = complex(zeros(n,m))
    if α > 0
        mm = N₋(m)
        A[:,1:mm] = reverse(CauchyPNO(n,mm,α,gridPts),dims=2)
    else
        mm = N₊(m)
        A[:,end-mm+1:end] = CauchyPNO(n,mm,α,gridPts)
    end
    return A
end

function *(C::CauchyOperator,domain::OscRational) #confused about how to do C+ without a BasisExpansion to call CauchyP...
    α = domain.α
    gd = domain.GD
    range = GridValues(gd)
    gridPts = gd.grid
    #working
    if C.o == 1.0
        if α == 0.
            return ConcreteOperator(domain,domain,BasicBandedOperator{ℤ,ℤ}(200,200, (i,j) -> CauchyConstantMatP(i,j))) #Need a better way to do this...
        else
            if α > 0
                Op1 = ConcreteOperator(domain,domain,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? complex(1.0) : 0.0im ))
            else
                Op1 = ConcreteOperator(domain,domain,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? complex(0.0) : 0.0im ))
            end
            Op2 = ConcreteOperator(domain,range,GenericEvaluationOperator{ℤ,𝔼}((n,m) -> BuildOperatorBlock(n,m,α,gridPts)))
            Op3 = Conversion(OscRational(gd,0.))
            return (Op1)⊘(Op3*Op2)
        end
    elseif C.o == -1.0
        if domain.α == 0.
            return ConcreteOperator(domain,domain,BasicBandedOperator{ℤ,ℤ}(200,200, (i,j) -> CauchyConstantMatM(i,j))) #Need a better way to do this...
        else
            # C+ - C- = I => C- = C+ - I (only affects oscillatory piece)
            if α > 0
                Op1 = ConcreteOperator(domain,domain,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? complex(0.0) : 0.0im )) #subtracted I from what is in C+ above
            else
                Op1 = ConcreteOperator(domain,domain,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? complex(-1.0) : 0.0im )) #subtract I from what is in C+ above
            end
            Op2 = ConcreteOperator(domain,range,GenericEvaluationOperator{ℤ,𝔼}((n,m) -> BuildOperatorBlock(n,m,α,gridPts))) #nonoscillatory piece is not affected so stays the same
            Op3 = Conversion(OscRational(gd,0.))
            return (Op1)⊘(Op3*Op2)
        end
    end
end

function Base.sum(f::BasisExpansion{OscRational}) #name change
    α = convert(Float64,f.basis.α)
    Laguerre = Lag(N₋(length(f.c)),2*abs(α))
    j_vals = -N₋(length(f.c)):N₊(length(f.c)) #determine j's
    sum = 0 #initialize integral
    for i=1:length(f.c)
         j = j_vals[i]
         #determine R_hat_j(α)
         if α == 0
             R_hat_j = -2*π*abs(j)
         elseif sign(j) == sign(α) #typo in AKNS (or at least discrepancy between AKNS and ref [26])
             R_hat_j = 0
         else
            R_hat_j = -4*π*Laguerre[abs(j) - 1 + 2] #exp(-abs(α)) already in Lag function so taking it out here
         end
         sum += f.c[i]*R_hat_j #add to integral value
    end
    return sum
end

function Base.conj(f::BasisExpansion{OscRational})
    α = -1*(f.basis.α)
    sp = OscRational(f.basis.GD,α)
    if isodd(length(f.c))
        coeffs = reverse(Base.conj(f.c))
    else
        coeffs = complex(zeros(length(f.c)+1))
        coeffs[2:end] = reverse(Base.conj(f.c))
    end
    return BasisExpansion(sp,coeffs)
end

function dot(f::BasisExpansion{OscRational},g::BasisExpansion{OscRational})
    return Base.sum(Multiplication(f)*Base.conj(g))
end