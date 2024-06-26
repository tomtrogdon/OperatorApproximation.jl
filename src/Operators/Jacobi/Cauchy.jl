function *(C::CauchyTransform,domain::Jacobi)
    a = domain.GD.D.a
    b = domain.GD.D.b
    α = domain.α
    β = domain.β
    gd = JacobiMappedInterval(a,b,α,β)
    range = Hardy(Exterior(gd)) # Just to ensure the right weight is encoded
    ConcreteOperator(domain,range,BasicBandedOperator{ℕ₊,ℕ₊}(0,0, (i,j) -> i == j && i >= 0 ? complex(1.0) : 0.0im ))
end