function d(j,λ)
    j*sqrt(2*(j + 2λ)/j*(λ+1)/(2λ+1))
end

function *(D::Derivative,domain::Ultraspherical)
    if D.order == 1
        range = Ultraspherical(domain.λ+1,domain.GD)
        dom = domain.GD.D
        ConcreteLazyOperator(domain,range,SingleBandedOperator(-1,1, (i,j) -> 2/(dom.b-dom.a)*d(j-1,domain.λ)))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end