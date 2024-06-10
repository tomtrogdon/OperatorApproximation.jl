function DerivMatrix(i,j)
    α = 2; #WE NEED TO FIND A WORKAROUND TO HARD CODING THIS!!!!!!!!!
    if (i == (j - 1)) || (i == (j + 1))
        return -1im*j/2
    elseif (i == j)
        return 1im*((j)+α)
    else 
        return 0im
    end
end

#This builds the derivative operator and applies it appropriately
function *(D::Derivative,domain::OscRational)
    if D.order == 1
        range = domain
        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(1,1,DerivMatrix))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end
