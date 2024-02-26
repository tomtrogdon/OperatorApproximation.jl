function *(M::Multiplication,sp::Ultraspherical)
    if typeof(M.f.basis) <: Ultraspherical && isconvertible(M.f.basis,sp)
        a, b = Jacobi_ab(M.f.basis.λ - 1/2,M.f.basis.λ- 1/2)
        α, β = Jacobi_ab(sp.λ - 1/2,sp.λ- 1/2)
        f = n -> OPMultiplication(a,b,α,β,M.f.c,sparse(I,n,n))[1:n,1:n] |> sparse
        m = length(M.f.c)
        Op = SemiLazyBandedOperator(m,m,f,f(5))
    else 
        1 + 1 #TODO: just evaluate and expand, need transform
    end
    ConcreteLazyOperator(sp,sp,Op)
end

function *(M::Multiplication,f::BasisExpansion{T}) where T<: Ultraspherical
    a, b = Jacobi_ab(M.f.basis.λ - 1/2,M.f.basis.λ- 1/2)
    α, β = Jacobi_ab(f.basis.λ - 1/2,f.basis.λ- 1/2)
    OPMultiplication(a,b,α,β,M.f.c,f.c)
end