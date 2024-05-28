function DerivMatrix(i,j,α)
    i = i-1
    j = j-1
    if j == 0
        if (i == 0)
            return 1im*α
        else
            return 0im
        end
    elseif j == 1
        if (i == 0) || (i == 2)
            return -3*1im*j
        elseif i == 1
            return 1im*α
        else
            return -2*1im*j
        end
    elseif j == 2
        if (i == 1) || (i == 3)
            return 0im
        elseif i == 2
            return (3*1im*j)+(1im*α)
        else
            return 1im*j
        end
    elseif (i == (j - 1)) || (i == (j + 1))
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
#         dom = domain.GD.D
        ConcreteOperator(domain,range,BasicMatrixOperator{ℤ,ℤ}(DerivMatrix(i,j,α)))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end
