function *(C::CauchyTransform,domain::FixedGridValues)
    gd = domain.GD
    range = Hardy(Exterior(gd)) # Just to ensure the right weight is encoded
    ConcreteLazyOperator(domain,range,BasicBandedOperator{ℕ₊,ℕ₊}(0,0, (i,j) -> i == j && i >= 0 ? complex(1.0) : 0.0im ))
end