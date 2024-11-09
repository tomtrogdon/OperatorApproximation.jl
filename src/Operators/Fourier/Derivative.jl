function *(D::Derivative,domain::Fourier)
    if D.order == 1
        range = domain
        dom = domain.GD.D
        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? 2im*pi*j/(dom.b-dom.a) : 0im ))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end

function *(D::FloquetDerivative,domain::Fourier)
    if D.order == 1
        range = domain
        dom = domain.GD.D
        ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? 1im*D.μ + 2im*pi*j/(dom.b-dom.a) : 0im ))
    else
        FloquetDerivative(D.order-1,D.μ)*(FloquetDerivative(1,D.μ)*domain)
    end
end