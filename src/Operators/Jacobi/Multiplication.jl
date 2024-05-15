function *(M::Multiplication,sp::Jacobi)
    if typeof(M.f) <: Function
        a = sp.GD.D.a
        b = sp.GD.D.b
        GD = JacobiMappedInterval(a,b,-0.5,-0.5)
        ff = BasisExpansion(M.f,Jacobi(-0.5,-0.5,GD)) |> chop
    else 
        ff = M.f
    end
    
    if typeof(ff.basis) <: Jacobi && ff.basis.GD.D == sp.GD.D
        a, b = Jacobi_ab(ff.basis.α,ff.basis.β)
        α, β = Jacobi_ab(sp.α,sp.β)
        f = n -> OPMultiplication(a,b,α,β,ff.c,sparse(I,n,n))[1:n,1:n] |> sparse
        m = length(ff.c)
        Op = SemiLazyBandedOperator{ℕ₊,ℕ₊}(m,m,f,f(5))
    else 
        1 + 1 #TODO: just evaluate and expand, need transform
    end
    ConcreteOperator(sp,sp,Op)
end