function rhdomain(endpoints::Matrix)
     temp = [(endpoints[i,1],endpoints[i,2]) for i=1:size(endpoints,1)]
     return ⊕([Legendre(t...) for t in temp]...)
end

function _rhrange(D::Basis)
    a = D.GD.D.a
    b = D.GD.D.b
    DirectedLobattoMappedInterval(a,b)
end

function mult2x2(A,B) # should be generalized
    a11 = z -> A[1,1](z)*B[1,1](z) + A[1,2](z)*B[2,1](z)
    a12 = z -> A[1,1](z)*B[1,2](z) + A[1,2](z)*B[2,2](z)
    a21 = z -> A[2,1](z)*B[1,1](z) + A[2,2](z)*B[2,1](z)
    a22 = z -> A[2,1](z)*B[1,2](z) + A[2,2](z)*B[2,2](z)
    [a11 a12; a21 a22]
end

function mvf2mof(f,n,m) # a bit lazy, tbh
    out = Matrix{Any}(nothing,n,m)
    for i = 1:n
        for j = 1:m
            out[i,j] = z -> f(z)[i,j]
        end
    end
    convert(Matrix{Function},out)
end

function mofeval(f,z)
    map(x -> x(z),f)
end

function rhrange(D::DirectSum)
    ⊕(GridValues.(_rhrange.(D.bases))...)
end

function rhrange(D::Basis)
    GridValues(_rhrange(D))
end

function RHmult(J::Matrix) 
    m = size(J,1) # m is size of RHP
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            Js[j,i] = Multiplication(J[i,j])
        end
    end
    convert(Matrix{Multiplication},Js)
end

# function RHmult(Js::Vector{T}) where T # J is a vector of scalar-valued functions

# end

function rhmult(J::Vector{T}) where T <: Matrix # J is a vector of matrices of scalar-valued functions
    if length(J) == 1
        return RHmult(J[1])
    end
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            g = [JJ[i,j] for JJ in J]
            Js[j,i] = BlockDiagonalAbstractOperator(Multiplication.(g))
        end
    end
    convert(Matrix{BlockDiagonalAbstractOperator},Js)
end

function rhrhs(J::Vector{T},c) where T <: Matrix # J is a vector of matrices of scalar-valued functions
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Vector{Any}(nothing,m)
    id = Matrix(I,m,m)
    for i = 1:m
        Js[i] = [ z -> (c*(ComplexF64.(map(x -> x(z),JJ)) - id))[i] for JJ in J]
    end
    convert(Vector{Vector},Js)
end

struct RHSolver
    S::ConcreteLazyOperator
    jumps
end

struct RHP
    Γ::Matrix
    J::Vector
end

contourplot(rhp::RHP;kwargs...) = contourplot(rhdomain(rhp.Γ);kwargs...)

function truncateRHP(Jsamp,J,Γ,tol,n)
    Gsamp = copy(Jsamp)
    G = copy(J)
    doms = Γ |> copy
    k = size(doms,1)
    
    doms = [doms[i,:] for i = 1:size(doms,1)]
    i = 0
    while i < k
        i += 1
        gd = LobattoMappedInterval(doms[i][1],doms[i][2])
        N = round(Int,n*arclength(gd)) 
        x = gd.D.map.(gd.grid(N))
        vals = abs.(Gsamp[i].(x))
        j = 1
        if vals[1] < tol
            for v in vals
                if v > tol
                    break
                end
                j += 1
            end
        end
        a = x[max(1,j-1)]
        l = length(vals)
        if vals[end] < tol
            for v in reverse(vals)
                if v > tol
                    break
                end
                l -= 1
            end
        end
        b = x[min(length(vals),l+1)]
        if j == length(vals) + 1 || l == 0
            deleteat!(doms,i)
            deleteat!(G,i)
            deleteat!(Gsamp,i)
            k -= 1
            i -= 1
        else
            doms[i] = [a, b]
        end
    end
    doms = [transpose(x) for x in doms]
    G, vcat(doms...)
end

function adapt(rhp::RHP,j,ϵ::Float64)
    J, Σ = truncateRHP(j,rhp.J,rhp.Γ,ϵ,100)
    RHP(Σ,J)
end

function (R::RHSolver)(c,n)
    b = vcat(RHrhs(R.jumps,c)...)
    u = \(R.S,b,n)
    k = length(R.jumps)
    m = length(c)
    if k == 1
        return [u[i] for i=1:m]
    end
    [u[(i-1)*k+1:i*k] for i=1:m]
end

function RHSolver(rhp::RHP)
    m = size(rhp.J[1],1) # size of RHP
    k = size(rhp.Γ,1) # number of intervals
    dom = rhdomain(rhp.Γ)
    ran = rhrange(dom)
    ℰ⁻ = BoundaryValue(-1,ran)
    ℰ⁺ = BoundaryValue(+1,ran)
    𝒞 = BlockAbstractOperator(CauchyTransform(),k,k)
    𝒞⁺ = ℰ⁺*𝒞
    𝒞⁻ = ℰ⁻*𝒞
    ℳ = RHmult(rhp.J)
    ℳ𝒞⁻ = matrix2BlockOperator(ℳ.*fill(𝒞⁻,m,m))
    𝒞⁺ = diagm(fill(𝒞⁺,m))
    dom = ⊕([dom for i = 1:m]...)
    ran = ⊕([ran for i = 1:m]...)
    S = (-ℳ𝒞⁻ + 𝒞⁺)*dom
    RHSolver(S,rhp.J)
end

### Vector "optimized" versions... that are slower... ###
struct RHSolverVec
    𝒞⁺::ConcreteLazyOperator
    𝒞⁻::ConcreteLazyOperator
    ℳ::ConcreteLazyOperator
    jumps
    range
    domain
end

function RHSolverVec(intervals::Matrix,jumps::Vector)
    m = size(jumps[1],1) # size of RHP
    k = size(intervals,1) # number of intervals
    dom = rhdomain(intervals)
    ran = rhrange(dom)
    ℰ⁻ = BoundaryValue(-1,ran)
    ℰ⁺ = BoundaryValue(+1,ran)
    𝒞 = BlockAbstractOperator(CauchyTransform(),k,k)
    𝒞⁺ = ℰ⁺*𝒞
    𝒞⁻ = ℰ⁻*𝒞
    ℳ = rhmult(jumps)
    ℳ = matrix2BlockOperator(map(x -> diagm(x.Ops),ℳ))
    #ℳ = diagm.(ℳ.Ops)
    RHSolverVec(𝒞⁺*dom,𝒞⁻*dom,ℳ*(ran ⊕ ran),jumps, ran, dom) 
end

# Only use for multiple contours
function (R::RHSolverVec)(c,n::Int64)
    ns1 = divide_DOF(R.range,n)
    m = length(c)
    k = length(ns1)
    ranges = vcat(fill(R.range.bases,m)...)
    domains = vcat(fill(R.domain.bases,m)...)
    ns = vcat(ns1,ns1)
    b = vcat(rhrhs(R.jumps,c)...)
    rhss = []
    for i = 1:length(ns)
        temp = BasisExpansion(b[i],ranges[i],ns[i])
        push!(rhss,temp.c)
    end
    b = vcat(rhss...)
    𝒞⁻ = Matrix(R.𝒞⁻,ns1,ns1) |> sparse
    𝒞⁺ = Matrix(R.𝒞⁺,ns1,ns1) |> sparse 
    𝒞⁻ = blockdiag(fill(𝒞⁻,m)...)
    𝒞⁺ = blockdiag(fill(𝒞⁺,m)...)
    ℳ = Matrix(R.ℳ,ns,ns)
    sol = (𝒞⁺ - ℳ*𝒞⁻)\b
    parted_sol = part_vec(sol,ns)
    u = ⊕(BasisExpansion.(domains,parted_sol)...)
    [u[(i-1)*k+1:i*k] for i=1:m]
end

function dilog(z)
    if abs(z) <= 3/4
        sum = 0.0
        Z = z
        for i = 1:95
            sum += Z/(i)^2
            Z *= z
        end
        return sum
    elseif abs(z) >= 4/3
        return -pi^2/6 - (log(-z |> complex))^2/2 - dilog(1/z)
    elseif abs(1-z) <= 3/4 || abs(1-z) >= 4/3
        return pi^2/6 - log(z)*log(1-z) - dilog(1-z)
    else
        w = sqrt(z)
        2*(dilog(w) + dilog(-w))
    end
end

function endpoint_list(dd)
    c = []
    for j = 1:length(dd)
        d = dd[j]
        an = d.GD.D.map(ArgNum(1.0,1.0,1.0*pi))
        push!(c,(an.z,an.θ,j,1))
        an = d.GD.D.map(ArgNum(-1.0,1.0,0.0))
        push!(c,(an.z,an.θ,j,-1))
    end
    c
end

function peel_endpoint(c)
    cc = [c[1]]
    inds = [1]
    ccpy = copy(c)
    for i = 2:length(c)
        cccc = c[i]
        if abs(c[1][1]-cccc[1]) < 1e-14
            push!(cc,cccc)
            push!(inds,i)
        end
    end
    deleteat!(ccpy,inds)
    (ccpy,cc)
end

function endpoint_check(ept,J)
    epts = sort(ept; lt = (x,y) -> x[2] < y[2])
    z = epts[1][1]
    σ = epts[1][4]
    At = mofeval(J[epts[1][3]],z)
    if σ == 1
        At = inv(At)
    end
    A = At
    for i = 2:length(epts)
        σ = epts[i][4]
        At = mofeval(J[epts[i][3]],z)
        if σ == 1
            At = inv(At)
        end
        A = A*At
    end
    (z, A)
end

function rhwellposed(rhp::RHP)
    el = rhp.Γ |> rhdomain |> endpoint_list
    out = []
    while length(el) > 0
        el, ept = peel_endpoint(el)
        push!(out,endpoint_check(ept,rhp.J))
    end
    out
end

function rhplot(rhp::RHP;kwargs...)
    # need to extend for larger RHPs
    dom = rhp.Γ |> rhdomain
    p0 = domainplot(dom;kwargs...)
    ran = rhrange(dom)
    N = 100
    plts = [p0]
    y = 0
    for i = 1:size(rhp.Γ,1)
        d = dom[i]
        x = d.GD.grid(N)
        z = d.GD.D.map.(x)
        y = vcat(map( x -> reshape(mofeval(rhp.J[i],x),1,:), z)...)
        p1 = plot(x,y[:,1] |> real;legend = false, kwargs...)
        plot!(p1,x,y[:,1] |> imag;legend = false, kwargs...)
        for k = 2:size(y,2)
            plot!(p1,x,y[:,k] |> real;legend = false, kwargs...)
            plot!(p1,x,y[:,k] |> imag;legend = false, kwargs...)
        end
        push!(plts,p1)
    end
    plts
end