function *(D::Derivative,domain::Laurent)
    if D.order == 1
        range = domain
        dom = domain.GD.D
        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(-1,1, (i,j) -> i + 1 == j ? j/dom.rad : 0im ))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end