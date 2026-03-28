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

function rhmult_jump(J::Matrix) 
    m = size(J,1) # m is size of RHP
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            Js[j,i] = Multiplication(J[i,j])
        end
    end
    convert(Matrix{Multiplication},Js)
end

# function rhmult(Js::Vector{T}) where T # J is a vector of scalar-valued functions

# end

function rhmult_jump(J::Vector{T}) where T <: Matrix # J is a vector of matrices of scalar-valued functions
    if length(J) == 1
        return rhmult_jump(J[1])
    end
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            g = [JJ[i,j] for JJ in J]
            Js[j,i] = BlockDiagonalAbstractOperator(Multiplication.(g)) # transpose is here
        end
    end
    convert(Matrix{BlockDiagonalAbstractOperator},Js)
end

function rhmult_res(J::Vector{T}) where T <: Matrix # J is a vector of matrices (res conds)
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            g = [JJ[i,j] for JJ in J]
            Js[j,i] = Multiplication(g) # transpose is here
        end
    end
    Js
end

function rhrhs_jump(J::Vector{T},c) where T <: Matrix # J is a vector of matrices of scalar-valued functions
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Vector{Any}(nothing,m)
    id = Matrix(I,m,m)
    for i = 1:m
        Js[i] = [ z -> (c*(ComplexF64.(map(x -> x(z),JJ)) - id))[i] for JJ in J]
    end
    convert(Vector{Vector},Js)
end

function rhrhs_res(J::Vector{T},c) where T <: Matrix # J is a vector of matrices (res conds)
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Vector{Any}(nothing,m)
    for i = 1:m
        Js[i] = [(c*JJ)[i] for JJ in J]
    end
    convert(Vector{Vector},Js)
end

struct RHSolver
    S::ConcreteOperator
    jumps
    res
end
RHSolver(S,jumps) = RHSolver(S,jumps,[])

struct RHP{T <: Union{Matrix,Vector}}
    Γ::T
    J::Vector
    P::Vector
    R::Vector
end
RHP(Γ,J) = RHP(Γ,J,[],[])

### Adaptive stuff ###
# TODO: Adapt for residues
function truncateRHP(Jsamp,J,Γ,P,R,tol,n,nvec=[])
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
        x = gd.D.map.(gd.grid(N+2))
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
            if !isempty(nvec)
                deleteat!(nvec,i)
            end
            k -= 1
            i -= 1
        else
            doms[i] = [a, b]
        end
    end
    doms = [transpose(x) for x in doms]

    poles = copy(P)
    res = copy(R)
    r = length(poles)
    for j = r:-1:1
        if norm(res[j]) < tol
            deleteat!(poles,j)
            deleteat!(res,j)
        end
    end

    if isempty(nvec)
        return G, vcat(doms...), poles, res
    end
    G, vcat(doms...), poles, res, nvec
end
#
function adapt(rhp::RHP,j,ϵ::Float64)
    J, Σ, P, R = truncateRHP(j,rhp.J,rhp.Γ,rhp.P,rhp.R,ϵ,100)
    RHP(Σ,J,P,R)
end

function adapt(rhp::RHP,j,ϵ::Float64,nvec::Vector{Int})
    J, Σ, P, R, nvec = truncateRHP(j,rhp.J,rhp.Γ,rhp.P,rhp.R,ϵ,100,nvec)
    RHP(Σ,J,P,R), nvec
end
########

function vcat_rhs(jumps,res,c)
    if length(res) == 0
        return vcat(rhrhs_jump(jumps,c)...)
    elseif length(jumps) == 0
        return rhrhs_res(res,c)
    else
        temp = rhrhs_jump(jumps,c)
        temp2 = rhrhs_res(res,c)
        return vcat(vcat.(map(x -> [x],temp2),temp)...)
    end
end

function (R::RHSolver)(c,n)
    b = vcat_rhs(R.jumps,R.res,c)
    u = \(R.S,b,n)
    k = length(R.jumps) + (length(R.res) > 0 ? 1 : 0)
    m = length(c)
    if k == 1
        return [u[i] for i=1:m]
    end
    return [u[(i-1)*k+1:i*k] for i=1:m]
end

function (R::RHSolver)(c::Tuple,n)
    b = map(c -> vcat_rhs(R.jumps,R.res,c), c)
    u = \(R.S,b,n)
    k = length(R.jumps) + (length(R.res) > 0 ? 1 : 0)
    m = length(c[1])
    q = length(c)
    out = []
    for j = 1:q
        if k == 1
            push!(out, [u[j][i] for i=1:m])
        else
            push!(out, [u[j][(i-1)*k+1:i*k] for i=1:m])
        end
    end
    return out
end

function (R::RHSolver)(c,n::Vector)
    if length(n) != length(R.jumps)
        error("Length of n must equal number of contours plus one if residues are present.")
    end
    b = vcat_rhs(R.jumps,R.res,c)
    m = length(c)
    if length(R.res) > 0
        nvec = repeat(vcat(length(R.res),n),m)
    else
        nvec = repeat(n,m)
    end
    u = \(R.S,b,nvec)
    k = length(R.jumps) + (length(R.res) > 0 ? 1 : 0)
    if k == 1
        return [u[i] for i=1:m]
    end
    return [u[(i-1)*k+1:i*k] for i=1:m]
end

function (R::RHSolver)(c::Tuple,n::Vector)
    if length(n) != length(R.jumps)
        error("Length of n must equal number of contours.")
    end
    b = map(c -> vcat_rhs(R.jumps,R.res,c), c)
    m = length(c[1])
    if length(R.res) > 0
        nvec = repeat(vcat(length(R.res),n),m)
    else
        nvec = repeat(n,m)
    end
    u = \(R.S,b,nvec)
    q = length(c)
    out = []
    k = length(R.jumps) + (length(R.res) > 0 ? 1 : 0)
    for j = 1:q
        if k == 1
            push!(out, [u[j][i] for i=1:m])
        else
            push!(out, [u[j][(i-1)*k+1:i*k] for i=1:m])
        end
    end
    return out
end

function RHSolver(rhp::RHP{T}) where T <: Vector
    m = size(rhp.R[1],1) # size of RHP
    dom = FixedGridValues(Grid(rhp.P))
    ran = dom
    ℰ⁻ = BoundaryValue(-1,ran)
    ℰ⁺ = Residue(dom)
    𝒞 = BlockAbstractOperator(CauchyTransform(),1,1)
    𝒞⁺ = ℰ⁺*𝒞
    𝒞⁻ = ℰ⁻*𝒞
    ℳ = rhmult_res(rhp.R)
    ℳ𝒞⁻ = matrix2BlockOperator(ℳ.*fill(𝒞⁻,m,m))
    𝒞⁺ = diagm(fill(𝒞⁺,m))
    dom = ⊕([dom for i = 1:m]...)
    ran = ⊕([ran for i = 1:m]...)
    S = (-ℳ𝒞⁻ + 𝒞⁺)*dom
    if length(rhp.P) > 0
        return RHSolver(S,rhp.J,rhp.R)
    end
    RHSolver(S,rhp.J)
end

## TODO:  Allow empty matrix of contours.
function RHSolver(rhp::RHP{T}) where T <: Matrix
    m = size(rhp.J[1],1) # size of RHP
    k = size(rhp.Γ,1) # number of intervals
    dom = rhdomain(rhp.Γ)
    ran = rhrange(dom)
    if length(rhp.P) > 0
        k += 1
        resdom = FixedGridValues(Grid(rhp.P))
        dom = resdom ⊕ dom
    end
    ℰ⁺ = BoundaryValue(+1,ran)
    if length(rhp.P) > 0
        ran = resdom ⊕ ran
    end
    ℰ⁻ = BoundaryValue(-1,ran)
    if length(rhp.P) > 0
        ℰ⁺ = Residue(resdom) ⊕ ℰ⁺
    end
    𝒞 = BlockAbstractOperator(CauchyTransform(),k,k)
    𝒞⁺ = ℰ⁺*𝒞
    𝒞⁻ = ℰ⁻*𝒞
    ℳ = rhmult_jump(rhp.J)
    if length(rhp.P) > 0
        ℳ = rhmult_res(rhp.R) .⊕ ℳ
    end
    ℳ𝒞⁻ = matrix2BlockOperator(ℳ.*fill(𝒞⁻,m,m))
    𝒞⁺ = diagm(fill(𝒞⁺,m))
    dom = ⊕([dom for i = 1:m]...)
    ran = ⊕([ran for i = 1:m]...)
    S = (-ℳ𝒞⁻ + 𝒞⁺)*dom
    if length(rhp.P) > 0
        return RHSolver(S,rhp.J,rhp.R)
    end
    RHSolver(S,rhp.J)
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

### For well-posed check ###
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
#
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
#
function endpoint_check(ept,J)
    epts = sort(ept; lt = (x,y) -> x[2] < y[2])
    z = epts[1][1]
    σ = epts[1][4]
    At = complex.(mofeval(J[epts[1][3]],z))
    if σ == 1
        At = inv(At)
    end
    A = At
    for i = 2:length(epts)
        σ = epts[i][4]
        At = complex.(mofeval(J[epts[i][3]],z))
        if σ == 1
            At = inv(At)
        end
        A = A*At
    end
    (z, A)
end
#
function rhwellposed(rhp::RHP)
    if size(rhp.Γ,1) == 1
        out = []
        z = rhp.Γ[1,1]
        y = mofeval(rhp.J[1],z)
        push!(out,(z,y))
        z = rhp.Γ[1,2]
        y = mofeval(rhp.J[1],z)
        push!(out,(z,y))
        return out
    end
    el = rhp.Γ |> rhdomain |> endpoint_list
    out = []
    while length(el) > 0
        el, ept = peel_endpoint(el)
        push!(out,endpoint_check(ept,rhp.J))
    end
    out
end

### Jacobi-based multi-interval RHP ###

struct JacobiRHP
    intervals::Matrix              # K×2: row k = [aₖ  bₖ]  (real or complex endpoints)
    Js::Vector                     # K m×m jump matrices
    αs::Union{Nothing,Matrix}
    βs::Union{Nothing,Matrix}      # K×m Jacobi exponents (nothing = auto-compute)
end

JacobiRHP(intervals, Js) = JacobiRHP(intervals, Js, nothing, nothing)

struct JacobiRHSolver
    S::ConcreteOperator
    Uis::Vector   # K inverse-eigenvector matrices
    K::Int
    m::Int
    Js::Vector    # K jump matrices
end

function JacobiRHSolver(rhp::JacobiRHP)
    intervals = rhp.intervals
    K = size(intervals, 1)

    # Normalise: wrap constant matrices as constant functions
    Js = [isa(Jk, Function) ? Jk : (z -> Jk) for Jk in rhp.Js]

    # Evaluate each jump at the interval midpoint for eigendecomposition
    mids   = [(intervals[k,1] + intervals[k,2]) / 2 for k = 1:K]
    Js_ref = [Js[k](mids[k]) for k = 1:K]
    m = size(Js_ref[1], 1)

    Es  = [eigen(Jk) for Jk in Js_ref]
    Uis = [inv(E.vectors) for E in Es]
    αs  = if rhp.αs !== nothing
        [rhp.αs[k, :] for k = 1:K]
    else
        [log.(E.values) ./ (2im*π) |> real for E in Es]
    end

    βs = if rhp.βs == nothing
        -αs
    else
        [rhp.βs[k, :] for k = 1:K]
    end

    # When αs[k][i] = 0... don't use Legendre
    function αeff(k, i)
        abs(αs[k][i]) > 1e-10 && return αs[k][i]
        for j = 1:m
            abs(αs[k][j]) > 1e-10 && return αs[k][j]
        end
        return 0.0
    end

    function βeff(k, i)
        abs(βs[k][i]) > 1e-10 && return βs[k][i]
        for j = 1:m
            abs(βs[k][j]) > 1e-10 && return βs[k][j]
        end
        return 0.0
    end

    gd(k, i) = JacobiMappedInterval(intervals[k,1], intervals[k,2], αeff(k,i), βeff(k,i))
    sp(k, i) = Jacobi(αeff(k,i), βeff(k,i), gd(k, i))
    gv(k)    = GridValues(gd(k, 1))

    dom   = ⊕([⊕([sp(k, i) for k = 1:K]...) for i = 1:m]...)
    ran_1 = ⊕([gv(k) for k = 1:K]...)

    ℰ⁺ = BoundaryValue(+1, ran_1)
    ℰ⁻ = BoundaryValue(-1, ran_1)
    𝒞  = BlockAbstractOperator(fill(CauchyTransform(), K, K))
    𝒞⁺ = ℰ⁺ * 𝒞
    𝒞⁻ = ℰ⁻ * 𝒞

    Gps = Matrix{Any}(nothing, m, m)
    Gms = Matrix{Any}(nothing, m, m)
    for i = 1:m, j = 1:m
        Gps[j,i] = BlockAbstractOperator([Multiplication(z -> Uis[l][i,j])           for k = 1:K, l = 1:K])
        Gms[j,i] = BlockAbstractOperator([Multiplication(z -> (Uis[l]*Js[k](z))[i,j]) for k = 1:K, l = 1:K])
    end
    Mp = convert(Matrix{BlockAbstractOperator}, Gps) |> BlockAbstractOperator
    Mm = convert(Matrix{BlockAbstractOperator}, Gms) |> BlockAbstractOperator

    ℳm𝒞⁻ = matrix2BlockOperator(Mm.Ops .⊙ fill(𝒞⁻, m, m))
    ℳp𝒞⁺ = matrix2BlockOperator(Mp.Ops .⊙ fill(𝒞⁺, m, m))
    Op = ℳp𝒞⁺ - ℳm𝒞⁻

    JacobiRHSolver(Op * dom, Uis, K, m, Js)
end

function (R::JacobiRHSolver)(c, N::Integer)
    Rhs = Vector{Any}(undef, R.m * R.K)
    for i = 1:R.m, k = 1:R.K
        Rhs[(i-1)*R.K + k] = x -> (c * ComplexF64.(R.Js[k](x) - I))[i]
    end
    u  = \(R.S, Rhs, N)
    Cu = reshape([CauchyTransform()*u[i] for i in 1:length(u)], R.K, R.m)
    u = reshape([u[i] for i in 1:length(u)], R.K, R.m)
    β = [moment(s,1) for s in u]*1im/(2pi)
    α = [moment(s,0) for s in u]*1im/(2pi)
    for i = 1:R.K
        β[i,:] = transpose(R.Uis[i])*β[i,:]
        α[i,:] = transpose(R.Uis[i])*α[i,:]
    end
    β = sum(β,dims=1)
    α = sum(α,dims=1)

    function Φ(z)
        data = reshape(Cu(z), R.K, R.m)
        for i = 1:R.K
            data[i,:] = transpose(R.Uis[i])*data[i,:]
        end
        c + sum(data, dims=1)
    end
    Φ, α, β, R.Uis, u
end

### Generalized Jacobi RHP: eigenvector matrices are functions of z ###
### THIS IS CURRENT NOT WORKING ### 

struct GeneralizedJacobiRHP
    intervals::Matrix              # K×2: row k = [aₖ  bₖ]  (real or complex endpoints)
    Js::Vector                     # K m×m jump matrices
    αs::Union{Nothing,Matrix}
    βs::Union{Nothing,Matrix}      # K×m Jacobi exponents (nothing = auto-compute)
end

GeneralizedJacobiRHP(intervals, Js) = GeneralizedJacobiRHP(intervals, Js, nothing, nothing)

struct GeneralizedJacobiRHSolver
    S::ConcreteOperator
    Uis::Vector   # K inverse-eigenvector matrices
    K::Int
    m::Int
    Js::Vector    # K jump matrices
end

function GeneralizedJacobiRHSolver(rhp::GeneralizedJacobiRHP)
    intervals = rhp.intervals
    K = size(intervals, 1)

    # Normalise: wrap constant matrices as constant functions
    Js = [isa(Jk, Function) ? Jk : (z -> Jk) for Jk in rhp.Js]

    # Evaluate each jump at the interval midpoint for eigendecomposition
    mids   = [(intervals[k,1] + intervals[k,2]) / 2 for k = 1:K]
    Js_ref = [Js[k](mids[k]) for k = 1:K]
    m = size(Js_ref[1], 1)

    Es  = [eigen(Jk) for Jk in Js_ref]
    Uis = [inv(E.vectors) for E in Es]
    αs  = if rhp.αs !== nothing
        [rhp.αs[k, :] for k = 1:K]
    else
        [log.(E.values) ./ (2im*π) |> real for E in Es]
    end

    βs = if rhp.βs == nothing
        -αs
    else
        [rhp.βs[k, :] for k = 1:K]
    end

    # When αs[k][i] = 0... don't use Legendre
    function αeff(k, i)
        abs(αs[k][i]) > 1e-10 && return αs[k][i]
        for j = 1:m
            abs(αs[k][j]) > 1e-10 && return αs[k][j]
        end
        return 0.0
    end

    function βeff(k, i)
        abs(βs[k][i]) > 1e-10 && return βs[k][i]
        for j = 1:m
            abs(βs[k][j]) > 1e-10 && return βs[k][j]
        end
        return 0.0
    end

    gd(k, i) = JacobiMappedInterval(intervals[k,1], intervals[k,2], αeff(k,i), βeff(k,i))
    sp(k, i) = Jacobi(αeff(k,i), βeff(k,i), gd(k, i))
    gv(k)    = GridValues(gd(k, 1))

    dom   = ⊕([⊕([sp(k, i) for k = 1:K]...) for i = 1:m]...)
    ran_1 = ⊕([gv(k) for k = 1:K]...)

    ℰ⁺ = BoundaryValue(+1, ran_1)
    ℰ⁻ = BoundaryValue(-1, ran_1)
    𝒞  = BlockAbstractOperator(fill(CauchyTransform(), K, K))
    𝒞⁺ = ℰ⁺ * 𝒞
    𝒞⁻ = ℰ⁻ * 𝒞

    # Gps = Matrix{Any}(nothing, m, m)
    # Gms = Matrix{Any}(nothing, m, m)
    # for i = 1:m, j = 1:m
    #     Gps[j,i] = BlockAbstractOperator([Multiplication(z -> Uis[l][i,j])       for k = 1:K, l = 1:K])
    #     Gms[j,i] = BlockAbstractOperator([Multiplication(z -> (Uis[l]*Js[k](z))[i,j]) for k = 1:K, l = 1:K])
    # end

    Gps = Matrix{Any}(nothing, m, m)
    Gms = Matrix{Any}(nothing, m, m)
    id = Matrix(I,m,m)
    for i = 1:m, j = 1:m
        Gps[j,i] = BlockAbstractOperator([Multiplication(z -> id[i,j])       for k = 1:K, l = 1:K])
        Gms[j,i] = BlockAbstractOperator([Multiplication(z -> (Js[k](z))[i,j]) for k = 1:K, l = 1:K])
    end

    Ups = Matrix{Any}(nothing, m, m)
    Ums = Matrix{Any}(nothing, m, m)
    for i = 1:m, j = 1:m
        Ups[j,i] = BlockAbstractOperator([Multiplication(z -> Uis[l][i,j]*(k == l ? 1.0 : 0.0)) for k = 1:K, l = 1:K])
        Ums[j,i] = BlockAbstractOperator([Multiplication(z -> Uis[l][i,j]*(k == l ? 1.0 : 0.0)) for k = 1:K, l = 1:K])
    end

    Mp = convert(Matrix{BlockAbstractOperator}, Gps) |> matrix2BlockOperator
    Mm = convert(Matrix{BlockAbstractOperator}, Gms) |> matrix2BlockOperator
    Up = convert(Matrix{BlockAbstractOperator}, Ups) |> matrix2BlockOperator
    Um = convert(Matrix{BlockAbstractOperator}, Ums) |> matrix2BlockOperator

    CCm = matrix2BlockOperator(fill(𝒞⁻, m, m))
    CCp = matrix2BlockOperator(fill(𝒞⁺, m, m))

    # ℳm𝒞⁻ = Mm * (Um ⊙ CCm)
    # ℳp𝒞⁺ = Mp * (Up ⊙ CCp)

    ℳm𝒞⁻ = (Mm ⊙ CCm) * Um
    ℳp𝒞⁺ = (Mp ⊙ CCp) * Up

    Op = ℳp𝒞⁺ - ℳm𝒞⁻

    GeneralizedJacobiRHSolver(Op * dom, Uis, K, m, Js)
end

function (R::GeneralizedJacobiRHSolver)(c, N::Integer)
    Rhs = Vector{Any}(undef, R.m * R.K)
    for i = 1:R.m, k = 1:R.K
        Rhs[(i-1)*R.K + k] = x -> (c * ComplexF64.(R.Js[k](x) - I))[i]
    end
    u  = \(R.S, Rhs, N)
    Cu = reshape([CauchyTransform()*u[i] for i in 1:length(u)], R.K, R.m)
    u = reshape([u[i] for i in 1:length(u)], R.K, R.m)
    β = [moment(s,1) for s in u]*1im/(2pi)
    α = [moment(s,0) for s in u]*1im/(2pi)
    for i = 1:R.K
        β[i,:] = transpose(R.Uis[i])*β[i,:]
        α[i,:] = transpose(R.Uis[i])*α[i,:]
    end
    β = sum(β,dims=1)
    α = sum(α,dims=1)

    function Φ(z)
        data = reshape(Cu(z), R.K, R.m)
        for i = 1:R.K
            data[i,:] = transpose(R.Uis[i])*data[i,:]
        end
        c + sum(data, dims=1)
    end
    Φ, α, β, R.Uis, u
end