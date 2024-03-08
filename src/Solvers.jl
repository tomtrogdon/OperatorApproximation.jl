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

function \(L::ConcreteLazyOperator{D,R,T},b::Vector,ns::Vector{Int64},ms::Vector{Int64}) where {D<:Basis,R<:Basis,T<:BlockLazyOperator}
    ranges = bases(L.range)
    domains = bases(L.domain)
    dimflag = dim.(ranges) .< Inf
    Op = Matrix(L,ns,ms)

    rhss = []
    for i = 1:length(ns)
        if dimflag[i] # functional just push vector
            push!(rhss,b[i])
        else
            temp = BasisExpansion(b[i],ranges[i],ns[i])
            push!(rhss,temp.c)
        end
    end
    sol = Op\vcat(rhss...)
    parted_sol = part_vec(sol,ms)
    BasisExpansion.(domains,parted_sol)
end

function \(L::ConcreteLazyOperator{D,R,T},b::Vector,N::Integer) where {D<:Basis,R<:Basis,T<:BlockLazyOperator}
    ns, ms = divide_DOF(L,N,N)
    \(L,b,ns,ms)
end

function \(L::ConcreteLazyOperator{D,R,T},b,N::Integer) where {D<:Basis,R<:Basis,T<:LazyOperator}
    Op = Matrix(L,N,N)
    rhs = BasisExpansion(b,L.range,N)
    BasisExpansion(L.domain,Op\rhs.c)
end

function testconv(f::BasisExpansion)  ## TODO:  May be incorrect for Bi-infinite things
    if dim(f.basis) < Inf
        return true
    else
        return norm(f.c[end-4:end]) < tol  
    end
end

function testconv(f::Vector{T}) where T <: BasisExpansion
   testconv.(f) |> prod
end

function \(L::ConcreteOperator,b)
    if !(typeof(N) <: Integer)
        n = 32
        sol = \(L,b,n)
        bool = testconv(sol)
        while !bool
            n *= 2
            sol = \(L,b,n)
            bool = testconv(sol)
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
    functions::Vector{BasisExpansion}
end

function eigen(Op::AbstractOperator,sp::Basis,N::Integer)
    eigen([Op],sp,N)
end

function eigen(L::Vector{T},sp::Basis,N::Integer) where T <: AbstractOperator
    Ops = [J*sp for J in L]
    eigen(Ops,N)
end

function eigen(L::Vector{T},N::Integer) where T <: ConcreteOperator
    k = 0
    Ops = []
    for i = 1:length(L)
        if isfinite(rank(L[i]))
            k += rank(L[i])
            push!(Ops,Matrix(L[i],N))
        end
    end
    #Op = vcat(Ops...)
    #display(Op)
    ## Assume only one non Functional Operator, for now...
    ## Adjust # of rows if not...
    for i = 1:length(L)
        if !isfinite(rank(L[i]))
            push!(Ops,Matrix(L[i],N-k,N))
        end
    end
    Op = vcat(Ops...)
    E = eigen(Op |> Matrix)
    funcs = [BasisExpansion(L[1].domain,E.vectors[:,i]) for i = 1:size(E.vectors)[2]]
    ContinuousEigen(E.values,funcs)
end