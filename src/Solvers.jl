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

function \(L::ConcreteOperator{D,R,T},b::Vector,ns::Vector{Int64}) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    return \(L,b,ns,ns)
end

function \(L::ConcreteOperator{D,R,T},b::Tuple,ns::Vector{Int64}) where {D<:Basis,R<:Basis,T<:BlockMatrixOperator}
    return \(L,b,ns,ns)
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

# ─────────────────────────────────────────────────────────────────────
#  GMRES
#
#  gmres(L, b, ip; tol, maxiter, restart)
#
#  L        – operator (AbstractOperator or ConcreteOperator).
#             Applied as L*x, so it must support multiplication by a
#             BasisExpansion (or return a BasisExpansion when applied
#             to the initial guess basis).
#  b        – right-hand side (BasisExpansion)
#  ip       – inner product: ip(u, v) -> scalar.  Must be sesquilinear
#             (conjugate-linear in the first argument).
#  tol      – relative residual tolerance (default 1e-12)
#  maxiter  – maximum number of outer iterations (default 100)
#  restart  – Krylov subspace size before restart (default maxiter)
# ─────────────────────────────────────────────────────────────────────

function gmres(L, b::BasisExpansion, ip::Function;
               tol     = 1e-12,
               maxiter = 100,
               restart = maxiter,
               x0      = nothing)

    _ip_norm(v) = sqrt(real(ip(v, v)))

    # Initial guess: zero expansion in the same basis as b
    x = x0 === nothing ? BasisExpansion(b.basis, zero(b.c)) : x0

    b_norm = _ip_norm(b)
    if b_norm == 0
        return x
    end

    for _outer in 1:cld(maxiter, restart)
        r = b - L * x
        r_norm = _ip_norm(r)

        if r_norm / b_norm < tol
            return x
        end

        m = min(restart, maxiter)

        # Arnoldi basis vectors and upper Hessenberg matrix
        Q = Vector{BasisExpansion}(undef, m + 1)
        H = zeros(ComplexF64, m + 1, m)

        Q[1] = (1 / r_norm) * r

        # Givens rotation storage
        cs = zeros(ComplexF64, m)
        sn = zeros(ComplexF64, m)
        e1 = zeros(ComplexF64, m + 1)
        e1[1] = r_norm

        j_final = m
        for j in 1:m
            w = L * Q[j]

            # Modified Gram-Schmidt
            for i in 1:j
                H[i, j] = ip(Q[i], w)
                w = w - H[i, j] * Q[i]
            end
            H[j+1, j] = _ip_norm(w)

            if abs(H[j+1, j]) < 1e-14 * b_norm
                j_final = j
                break
            end

            Q[j+1] = (1 / H[j+1, j]) * w

            # Apply previous Givens rotations to new column
            for i in 1:j-1
                tmp        =  cs[i] * H[i, j] + sn[i] * H[i+1, j]
                H[i+1, j]  = -conj(sn[i]) * H[i, j] + cs[i] * H[i+1, j]
                H[i, j]    =  tmp
            end

            # Compute new Givens rotation
            denom = sqrt(abs2(H[j, j]) + abs2(H[j+1, j]))
            cs[j] = H[j, j]   / denom
            sn[j] = H[j+1, j] / denom

            H[j, j]   =  cs[j] * H[j, j] + sn[j] * H[j+1, j]
            H[j+1, j] = 0

            e1[j+1] = -conj(sn[j]) * e1[j]
            e1[j]   =  cs[j] * e1[j]

            if abs(e1[j+1]) / b_norm < tol
                j_final = j
                break
            end
        end

        # Solve upper-triangular system H[1:j,1:j] y = e1[1:j]
        j = j_final
        y = H[1:j, 1:j] \ e1[1:j]

        # Update solution
        for i in 1:j
            x = x + y[i] * Q[i]
        end

        # Check convergence after restart
        r = b - L * x
        if _ip_norm(r) / b_norm < tol
            return x
        end
    end

    return x
end
