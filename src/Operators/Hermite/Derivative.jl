function herm_d(i,j)
    if j == i + 1
        return -sqrt(i)
    elseif j == i - 1
        return sqrt(j)
    else
        return 0.0
    end
end

function *(D::Derivative,domain::HermiteFun)
    if D.order == 1
        ConcreteLazyOperator(domain,domain,BasicBandedOperator{ℕ₊,ℕ₊}(1,1,herm_d))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end