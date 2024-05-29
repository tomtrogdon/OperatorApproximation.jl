function DerivMatrix(i,j,α)
    if (i == (j - 1)) || (i == (j + 1))
        return -1im*j
    elseif (i == j)
        return 1im*((2*j)+α)
    else 
        return 0im
    end
end

#This builds the derivative operator and applies it appropriately
function *(D::Derivative,domain::Rational)
    if D.order == 1
        range = domain
        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(1,1,DerivMatrix(i,j,α)))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end
