function DerivMatrix(i,j,α)
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
    α = domain.α
    if D.order == 1
        range = domain
        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(1,1,(i,j) -> DerivMatrix(i,j,α)))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end

function WeightedDerivative(domain::OscRational) # really only makes sense for α = 0
    ConcreteOperator(domain,domain,BasicBandedOperator{ℤ,ℤ}(0,1, (i,j) -> i + 1 == j ? complex(j*1.0) : 0.0im ))
end