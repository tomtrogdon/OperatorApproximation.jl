struct Jacobi <: Basis
    α::Number
    β::Number
    GD::GridInterval
end

function Legendre(a,b)
    GD = JacobiMappedInterval(a,b,0.0,0.0)
    Jacobi(0.0,0.0,GD)
end

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::Jacobi) = ℕ₊

function dim(sp::Jacobi)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Jacobi
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: Jacobi
    BasisExpansion(f.basis,chop(f.c))
end
####################################
#####  Important to implement  #####
####################################
function sum(f::BasisExpansion{T}) where T <: Jacobi
    (f.basis.GD.D.b - f.basis.GD.D.a)*f.c[1]
end

function (P::BasisExpansion{Jacobi})(X::Number) # Clenshaw's algorithm
    n = P.c |> length
    α = P.basis.α
    β = P.basis.β
    x = P.basis.GD.D.imap(X)
    a,b = Jacobi_ab(α,β)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
end

function Jacobi_ab(a,b) #TODO: simplify evaluation
    bfun = n -> (a+b==-1 && n==0) ? √(2*a*b) :
        2*sqrt(n+1)*sqrt(n+a+1)*sqrt(n+b+1)*sqrt(n+a+b+1)/
        ((2n + a +b + 2)*sqrt(2n + a +b + 3)*sqrt(2n + a +b + 1))
    afun = n -> ((a+b==0 || a+b==-1) && n==0) ? (b-a)/(a+b+2) :
        (b^2 - a^2)/((2n + a +b + 2)*(2n + a +b))
    return (n -> n < 0 ? a : afun(n),n -> n < 0 ? b : bfun(n))  # this is not needed but lets us get access
                                                                # to the Jacobi parameters after the fact
end

function jacobi(a,b,n) # creates (n + 1) x (n+1) Jacobi matrix
   SymTridiagonal([a(i) for i in 0:n],[b(i) for i in 0:n-1])
end

function Interp_transform(a,b,n)
    E = jacobi(a,b,n) |> eigen
    return E.values, E.vectors*(Diagonal(E.vectors[1,:]))
end

function OPMultiplication(a::Function,b::Function,α::Function,β::Function,c::Vector,u::Union{Array,SparseMatrixCSC})
    ## α, β give the recurrence coefficients for the basis used
    # on the domain of the operator.  a, b give the recurrence coefficitions
    # for the basis used for the function doing the multiplication.       
    n = size(u)[1] + length(c) + 3
    J = jacobi(α,β,n-1)
    shape = size(u)
    shape = tuple(length(c)+ 3,shape[2:end]...)
    U = vcat(u,zeros(shape)) |> sparse # TODO: Better way to do this?
    p = U # p_0
    q = c[1]*p
    polder = p
    p = J*p - a(0)*p # compute p_1
    p /= b(0)
    q += c[2]*p
    pold = p
    for j = 3:length(c) # compute p_n
        p = J*pold - a(j-2)*pold - b(j-3)*polder
        p /= b(j-2)
        q += c[j]*p
        polder = pold
        pold = p
    end
    q |> sparse
end

function Gauss_quad(a,b,n)
    E = jacobi(a,b,n) |> eigen
    return E.values, abs2.(E.vectors[1,:])
end

function Interp_transform(f::Function,a,b,n)
    E = jacobi(a,b,n) |> eigen
    E.vectors*(Diagonal(E.vectors[1,:])*map(f,E.values))
end

function poly(a,b,n,x) # a naive use of the three-term recurrence
    if n == 1
        return fill(1.0,n)
    end
    p = fill(0.0,n)
    p[1] = 1.0 # p_0
    p[2] = x.*p[1] - a(0)*p[1] # compute p_1
    p[2] /= b(0)
    for j = 1:n-2 # compute p_n
        p[j+2] = x.*p[j+1] - a(j)*p[j+1] - b(j-1)*p[j]
        p[j+2] /= b(j)
    end
    p
end

function poly(a,b,n,z::Vector)
    vcat(map(zz -> poly(a,b,n,zz) |> transpose , z)...)
end

function JacobiW(a,b,x)
    c = (_₂F₁(1, -a, 2+b, -1))/(1+b) + (_₂F₁(1, -b, 2+a, -1))/(1+a)
    return 1/c*(1-x)^a*(1+x)^b
end

function JacobiWconst(a,b)
    c = (_₂F₁(1, -a, 2+b, -1))/(1+b) + (_₂F₁(1, -b, 2+a, -1))/(1+a)
    return 1/c
end

## from SingularIntegralEquation.jl
normalization(n::Int,α::Real,β::Real) = 2^(α+β)*gamma(n+α+1)*gamma(n+β+1)/gamma(2n+α+β+2)
#stieltjesjacobimoment(α::Real,β::Real,n::Int,z) =
    #(x = 2/(1-z);normalization(n,α,β)*HypergeometricFunctions.mxa_₂F₁(n+1,n+α+1,2n+α+β+2,x))
stieltjesjacobimoment(α::Real,β::Real,n::Int,z) =
    (x = 2/(1-z);HypergeometricFunctions.mxa_₂F₁(n+1,n+α+1,2n+α+β+2,x))/2
stieltjesjacobimoment(α::Real,β::Real,z) = stieltjesjacobimoment(α,β,0,z)

function matanh(z)
    1/2*(log(1+z) - log(-1+z))
end

function matanh_p(z) # limit from above for z > 1
    1/2*(log(1+0im+z) - log(-1+0im+z)) + 1im*pi
end

function matanh_m(z) # limit from below for z > 1
    1/2*(log(1+0im+z) - log(-1+0im+z))
end

function legendrestieltjes(z)
    return 1im/(4*pi)*(log(-complex(1)-z)-log(complex(1)-z))
end

function legendrestieltjes_pos(z)
    return 1im/(4*pi)*(log(1+z) - 1im*pi -log(1-z))
end

function legendrestieltjes_neg(z)
    return 1im/(4*pi)*(log(1+z) + 1im*pi -log(1-z))
end

function JacobiSeed(α,β)
    if α == 0.0 && β == 0.0
        return z -> legendrestieltjes(z)
    elseif α == 0.5 && β == 0.0 # log singularity at z = -1
        out = JacobiSeed(β,α)
        return z -> -out(-z)
    elseif α == 0.0 && β == 0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ?  3im/(8π)*(-2 + log(4)) + 3/(4*sqrt(2))*2*sqrt(2)*legendrestieltjes(z) :
         1/(8im*pi)*(6 - 3*sqrt(2)*sqrt(1+z + 0im)*matanh(sqrt(1+z+ 0im)/sqrt(2)))
    elseif α == -0.5 && β == 0.0
        out = JacobiSeed(β,α)
        z -> -out(-z)
    elseif α == 0.0 && β == -0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ? 1im*log(2.0)/(4π) + 0.5*legendrestieltjes(z)  :
        -1/(4im*sqrt(2)*pi)*1/sqrt(1+z +0im)*log((z-1)/(3 + complex(z) - 2*sqrt(2)*sqrt(1 + z |> complex)))
    else
        return z -> 1im/(2*pi)*stieltjesjacobimoment(α,β,z)
    end
end

function JacobiSeedPos(α,β)
    if α == 0.0 && β == 0.0
        return z -> legendrestieltjes_pos(z)
    elseif α == 0.5 && β == 0.0 # log singularity at z = -1
        out = JacobiSeedNeg(β,α)
        return z -> -out(-z)
    elseif α == 0.0 && β == 0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ?  3im/(8π)*(-2 + log(4)) + 3/(4*sqrt(2))*2*sqrt(2)*legendrestieltjes_pos(z) :
         1/(8im*pi)*(6 - 3*sqrt(2)*sqrt(1+z + 0im)*matanh_m(sqrt(1+z+ 0im)/sqrt(2)))
    elseif α == -0.5 && β == 0.0
        out = JacobiSeedNeg(β,α)
        return z -> -out(-z)
    elseif α == 0.0 && β == -0.5
        return z -> abs(z - 1) < 1e-14 ? 1im*log(2.0)/(4π) + 0.5*legendrestieltjes_pos(z)  :
        -1/(4im*sqrt(2)*pi)*1/sqrt(1+z)*(log(-(z-1)/(3 + z - 2*sqrt(2)*sqrt(1 + z))) - 1im*pi )
    else
        return z -> 1im/(2*pi)*stieltjesjacobimoment(α,β,z + 1im*eps())
    end
end

function JacobiSeedNeg(α,β)
    if α == 0.0 && β == 0.0
        return z -> legendrestieltjes_neg(z)
    elseif α == 0.5 && β == 0.0 # log singularity at z = -1
        out = JacobiSeedPos(β,α)
        z -> -out(-z)
    elseif α == 0.0 && β == 0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ?  3im/(8π)*(-2 + log(4)) + 3/(4*sqrt(2))*2*sqrt(2)*legendrestieltjes_neg(z) :
        1/(8im*pi)*(6 - 3*sqrt(2)*sqrt(1+z)*(matanh_m(sqrt(1+z)/sqrt(2)) + 1im*pi))
    elseif α == -0.5 && β == 0.0
        out = JacobiSeedPos(β,α)
        z -> -out(-z)
    elseif α == 0.0 && β == -0.5
        return z -> abs(z - 1) < 1e-14 ? 1im*log(2.0)/(4π) + 0.5*legendrestieltjes_neg(z)  :
        -1/(4im*sqrt(2)*pi)*1/sqrt(1+z)*(log(-(z-1)/(3 + z - 2*sqrt(2)*sqrt(1 + z))) + 1im*pi )
    else
        return z -> 1im/(2*pi)*stieltjesjacobimoment(α,β,z - 1im*eps())
    end
end