#This builds the derivative operator and applies it appropriately
function *(D::Derivative,domain::Rational)
    if D.order == 1
        range = domain
        dom = domain.GD.D
        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? 2im*pi*j/(dom.b-dom.a) : 0im ))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end
