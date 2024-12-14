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
#
function adapt(rhp::RHP,j,ϵ::Float64)
    J, Σ = truncateRHP(j,rhp.J,rhp.Γ,ϵ,100)
    RHP(Σ,J)
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

function RHSolver(rhp::RHP{T}) where T <: Vector
    m = size(rhp.R[1],1) # size of RHP
    # k = size(rhp.Γ,1) # number of intervals
    # dom = rhdomain(rhp.Γ)
    # ran = rhrange(dom)
    dom = FixedGridValues(Grid(rhp.P))
    ran = dom

    # if length(rhp.P) > 0
    #     k += 1
    #     resdom = FixedGridValues(Grid(rhp.P))
    #     dom = resdom ⊕ dom
    # end
    ℰ⁻ = BoundaryValue(-1,ran)
    #ℰ⁺ = BoundaryValue(+1,ran)
    ℰ⁺ = Residue(dom)
    # if length(rhp.P) > 0
    #     ran = resdom ⊕ ran
    # end
    # ℰ⁻ = BoundaryValue(-1,ran)
    # if length(rhp.P) > 0
    #     ℰ⁺ = Residue(resdom) ⊕ ℰ⁺
    # end
    𝒞 = BlockAbstractOperator(CauchyTransform(),1,1)
    𝒞⁺ = ℰ⁺*𝒞
    𝒞⁻ = ℰ⁻*𝒞
    ℳ = rhmult_res(rhp.R)
    # ℳ = rhmult_jump(rhp.J)
    # if length(rhp.P) > 0
    #     ℳ = rhmult_res(rhp.R) .⊕ ℳ
    # end
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
#########