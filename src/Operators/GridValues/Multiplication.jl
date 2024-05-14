function *(M::Multiplication,sp::GridValues)
    Op = GridMultiplication(M.f, n -> sp.GD.D.map(sp.GD.grid(n)))
    ConcreteOperator(sp,sp,Op)
end

function *(M::Multiplication,sp::FixedGridValues)
    if typeof(M.f) <: Function || typeof(M.f) <: BasisExpansion
        Op = GridMultiplication(M.f.(sp.pts),sp.pts)
        return ConcreteOperator(sp,sp,Op)
    else
        Op = GridMultiplication(M.f,sp.pts)
        return ConcreteOperator(sp,sp,Op)
    end
end

function Matrix(GM::GridMultiplication,n,m)
    if typeof(GM.f) <: Function
        if n == m
            return Diagonal(GM.f.(GM.grid(n))) |> sparse
        else
            return Diagonal(GM.f.(GM.grid(max(n,m))))[1:n,1:m] |> sparse
        end
    else
        A = Diagonal(GM.f)
        k = length(GM.f)
        if n == m && n == k
            return A
        elseif n <= k && m <= k
            return A[1:n,1:m] |> sparse
        else
            @warn "One dimension too large.  Returning smaller matrix."
            return A[1:min(n,k),1:min(m,k)] |> sparse
        end
    end
end
    