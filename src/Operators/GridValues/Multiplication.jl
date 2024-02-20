struct GridMultiplication <: DenseOperator # even though it is sparse...
    # it is simpler to treat grid multiplication as dense
    f::Function
    grid::Function
end

struct FixedGridMultiplication <: DenseOperator
    fvals::Vector
end

function *(M::Multiplication,sp::GridValues)
    Op = GridMultiplication(M.f, n -> sp.GD.D.map(sp.GD.grid(n)))
    ConcreteLazyOperator(sp,sp,Op)
end

function *(M::Multiplication,sp::FixedGridValues)
    Op = FixedGridMultiplication(M.f(sp.GD.D.map(sp.pts)))
    ConcreteLazyOperator(sp,sp,Op)
end

function Matrix(GM::GridMultiplication,n,m)
    if n == m
        return Diagonal(GM.f(GM.grid(n)))
    else
        return Diagonal(GM.f(GM.grid(max(n,m))))[1:n,1:m] |> sparse
    end
end

function Matrix(GM::FixedGridMultiplication,n,m)
    A = Diagonal(GM.fvals)
    k = length(GM.fvals)
    if n == m && n == k
        return A
    elseif n <= k && m <= k
        return A[1:n,1:m] |> sparse
    else
        @warn "One dimension too large.  Returning smaller matrix."
        return A[1:min(n,k),1:min(m,k)] |> sparse
    end
end
    