function *(D::Derivative,domain::Erf)
    if D.order == 1
        range = HermiteFun(domain.GD)
        ConcreteLazyOperator(domain,range,BasicBandedOperator{ℕ₊,ℕ₊}(0,0, (i,j) -> i == 1 && j == 1 ? 1.0 : 0.0))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end