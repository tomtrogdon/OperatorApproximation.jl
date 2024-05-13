function herm_d(i,j)
    # if i == 1 && j == 2
    #     return 0.5
    # elseif j == 1 && i == 2
    #     return 0.5
    # else
    if j == i + 1
        return sqrt(i)/2
    elseif j == i - 1
        return -sqrt(j)/2
    else
        return 0.0
    end
end

function *(D::Derivative,domain::HermiteFun)
    if D.order == 1
        ConcreteOperator(domain,domain,BasicBandedOperator{ℕ₊,ℕ₊}(1,1,herm_d))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end