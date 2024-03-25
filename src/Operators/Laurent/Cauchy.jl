## Both CauchyTransform and CauchyOperator implemented here

function *(C::CauchyTransform,domain::Laurent)
    range = Hardy(Interior(domain.GD))
    Cp = ConcreteLazyOperator(domain,range,BasicBandedOperator{ℤ,ℕ₊}(0,0, (i,j) -> i == j && i >= 0 ? complex(1.0) : 0.0im ))
    range = Hardy(Exterior(domain.GD))
    Cm = ConcreteLazyOperator(domain,range,BasicBandedOperator{ℤ,ℕ₋}(0,0, (i,j) -> i == j && i < 0 ? complex(-1.0) : 0.0im ))
    Cp ⊘ Cm
end