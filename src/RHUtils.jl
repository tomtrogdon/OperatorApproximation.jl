function RHdomain(endpoints::Matrix)
     temp = [(endpoints[i,1],endpoints[i,2]) for i=1:size(endpoints,1)]
     return ⊕([Legendre(t...) for t in temp]...)
end

function _rhrange(D::Basis)
    a = D.GD.D.a
    b = D.GD.D.b
    DirectedLobattoMappedInterval(a,b)
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

function RHrange(D::DirectSum)
    ⊕(GridValues.(_rhrange.(D.bases))...)
end

function RHrange(D::Basis)
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

function RHmult(J::Vector{T}) where T <: Matrix # J is a vector of matrices of scalar-valued functions
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

function RHrhs(J::Vector{T},c) where T <: Matrix # J is a vector of matrices of scalar-valued functions
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

function RHSolver(intervals::Matrix,jumps::Vector)
    m = size(jumps[1],1) # size of RHP
    k = size(intervals,1) # number of intervals
    dom = RHdomain(intervals)
    ran = RHrange(dom)
    ℰ⁻ = BoundaryValue(-1,ran)
    ℰ⁺ = BoundaryValue(+1,ran)
    𝒞 = BlockAbstractOperator(CauchyTransform(),k,k)
    𝒞⁺ = ℰ⁺*𝒞
    𝒞⁻ = ℰ⁻*𝒞
    ℳ = RHmult(jumps)
    ℳ𝒞⁻ = matrix2BlockOperator(ℳ.*fill(𝒞⁻,m,m))
    𝒞⁺ = diagm(fill(𝒞⁺,m))
    dom = ⊕([dom for i = 1:m]...)
    ran = ⊕([ran for i = 1:m]...)
    S = (-ℳ𝒞⁻ + 𝒞⁺)*dom
    RHSolver(S,jumps)
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
    dom = RHdomain(intervals)
    ran = RHrange(dom)
    ℰ⁻ = BoundaryValue(-1,ran)
    ℰ⁺ = BoundaryValue(+1,ran)
    𝒞 = BlockAbstractOperator(CauchyTransform(),k,k)
    𝒞⁺ = ℰ⁺*𝒞
    𝒞⁻ = ℰ⁻*𝒞
    ℳ = RHmult(jumps)
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
    b = vcat(RHrhs(R.jumps,c)...)
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