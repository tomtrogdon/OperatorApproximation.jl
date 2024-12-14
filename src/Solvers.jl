function part_vec(x,ns)
    if length(x) != sum(ns)
        @error "Length issue in partitioning"
        return 0
    end
    xcs = cumsum(ns)
    out = [x[1:xcs[1]]]
    for i = 2:length(xcs)
        push!(out,x[xcs[i-1]+1:xcs[i]])
    end
    out
end


function _rhs_vec_gen(ns,dimflag,b,ranges)
    rhss = []
    for i = 1:length(ns)
        if dimflag[i] # functional just push vector
            push!(rhss,b[i])
        else
            temp = BasisExpansion(b[i],ranges[i],ns[i])
            push!(rhss,temp.c)
        end
    end
    vcat(rhss...)
end

function \(L::ConcreteOperator{D,R,T},b::Vector,ns::Vector{Int64},ms::Vector{Int64}) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    ranges = bases(L.range)
    domains = bases(L.domain)
    dimflag = dim.(ranges) .< Inf
    Op = Matrix(L,ns,ms)
    # rhss = []
    # for i = 1:length(ns)
    #     if dimflag[i] # functional just push vector
    #         push!(rhss,b[i])
    #     else
    #         temp = BasisExpansion(b[i],ranges[i],ns[i])
    #         push!(rhss,temp.c)
    #     end
    # end
    rhss = _rhs_vec_gen(ns, dimflag, b, ranges)   
    sol = Op\rhss
    #sol = lu!(Op)\rhss
    parted_sol = part_vec(sol,ms)
    ⊕(BasisExpansion.(domains,parted_sol)...)
end

function \(L::ConcreteOperator{D,R,T},b::Tuple,ns::Vector{Int64},ms::Vector{Int64}) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    ranges = bases(L.range)
    domains = bases(L.domain)
    dimflag = dim.(ranges) .< Inf
    Op = Matrix(L,ns,ms)
    rhss = map(b -> _rhs_vec_gen(ns, dimflag, b, ranges), b)
    rhss = hcat(rhss...)
    sol = Op\rhss
    #sol = lu!(Op)\rhss
    out = []
    for i = 1:length(b)
        parted_sol = part_vec(sol[:,i],ms)
        u = ⊕(BasisExpansion.(domains,parted_sol)...)
        push!(out,u)
    end
    out
end

function \(L::ConcreteOperator{D,R,T},b::Vector,N::Integer) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    ns, ms = divide_DOF(L,N,N)
    \(L,b,ns,ms)
end

function \(L::ConcreteOperator{D,R,T},b::Tuple,N::Integer) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    ns, ms = divide_DOF(L,N,N)
    \(L,b,ns,ms)
end

function \(L::ConcreteOperator{D,R,T},b,N::Integer) where {D<:Basis,R<:Basis,T<:MatrixOperator}
    Op = Matrix(L,N,N)
    rhs = BasisExpansion(b,L.range,N)
    BasisExpansion(L.domain,Op\rhs.c)
end

function \(L::ConcreteOperator{D,R,T},b::Tuple,N::Integer) where {D<:Basis,R<:Basis,T<:MatrixOperator}
    Op = Matrix(L,N,N)
    rhss = map(b -> BasisExpansion(b,L.range,N).c,b)
    rhss = hcat(rhss...)
    #sol = Op\rhss
    sol = lu!(Op)\rhss
    out = []
    for i = 1:length(b)
        u = BasisExpansion(L.domain,sol[:,i])
        push!(out,u)
    end
    out
end

function testconv(f::Vector{T}) where T <: Number
   k =  min(4,length(f) ÷ 2)
   maximum(abs.(f[end-k:end])) < 1e-13
end

function \(L::ConcreteOperator,b)
    if !(typeof(N) <: Integer)
        n = 32
        sol = \(L,b,n)
        if typeof(sol) <: Vector
            bool = testconv.(sol) |> prod
        else
            bool = testconv(sol)
        end
        while !bool
            n *= 2
            sol = \(L,b,n)
            if typeof(sol) <: Vector
                bool = testconv.(sol) |> prod
            else
                bool = testconv(sol)
            end
        end
        return sol
    else
        return \(L,b,N)
    end
end

function \(L::AbstractOperator,b)
    (L*basis)\b
end

function \(L::AbstractOperator,b,basis::Basis,N::Integer)
    \(L*basis,b,N)
end

function \(L::AbstractOperator,b,N::Integer)
    \(L*basis,b,N)
end

struct ContinuousEigen
    values::Vector{Union{Float64,ComplexF64}}
    functions::Vector
end

## non-generalized, block problem
function eigen(L::ConcreteOperator{D,R,T},ns::Vector{Int64},ms::Vector{Int64}) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    E =  eigen(Matrix(L,ns,ms) |> Matrix)
    domains = bases(L.domain)
    vs = [part_vec(E.vectors[:,i],ns) for i in 1:size(E.vectors,2)]
    vs = [BasisExpansion.(domains, vsp) for vsp in vs]
    ContinuousEigen(E.values,vs)
end
#
function eigen(L::ConcreteOperator{D,R,T},N::Int64) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    ns, ms = divide_DOF(L,N,N)
    eigen(L,ns,ms)
end

## non-generalized, single problem
function eigen(L::ConcreteOperator{D,R,T},N::Int64) where {D<:Basis,R<:Basis,T}
    E =  eigen(Matrix(L,N,N) |> Matrix)
    vs = [BasisExpansion(L.domain, E.vectors[:,i]) for i in 1:size(E.vectors,2)]
    ContinuousEigen(E.values,vs)
end

function makeinf(x)
    if abs(x) > 1e12
        Inf
    else
        x
    end
end

## generalized, block problem
function eigen(L::ConcreteOperator{D,R,T},M::ConcreteOperator{D,R,T},ns::Vector{Int64},ms::Vector{Int64}) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    E =  eigen(Matrix(L,ns,ms) |> Matrix, Matrix(M,ns,ms) |> Matrix)
    domains = bases(L.domain)
    vs = [part_vec(E.vectors[:,i],ms) for i in 1:size(E.vectors,2)]
    vs = [BasisExpansion.(domains, vsp) for vsp in vs]
    ContinuousEigen(makeinf.(E.values),vs)
end
#
function eigen(L::ConcreteOperator{D,R,T},M::ConcreteOperator{D,R,T},N::Integer) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    ns, ms = divide_DOF(L,N,N)
    eigen(L,M,ns,ms)
end

## generalized, single problem
function eigen(L::ConcreteOperator{D,R,T},M::ConcreteOperator{D,R,T},N::Integer) where {D<:Basis,R<:Basis,T}
    E =  eigen(Matrix(L,N,N) |> Matrix, Matrix(M,N,N) |> Matrix)
    domains = bases(L.domain)
    p = x -> part_vec(x,ns)
    vs = [BasisExpansion.(L.domain, E.vectors[:,i]) for i in 1:size(E.vectors,2)]
    ContinuousEigen(makeinf.(E.values),vs)
end
