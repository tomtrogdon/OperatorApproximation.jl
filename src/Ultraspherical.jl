struct Ultraspherical <: Basis
    λ::Number
    GD::GridDomain
end

function iscompatible(US1::Ultraspherical,US2::Ultraspherical)
    US1.λ == US2.λ && US1.GD == US2.GD
end

function d(j,λ)
    j*sqrt(2*(j + 2λ)/j*(λ+1)/(2λ+1))
end

function *(D::Derivative,domain::Ultraspherical)
    if D.order == 1
        range = Ultraspherical(domain.λ+1,domain.GD)
        dom = domain.GD.D
        ConcreteOperator(domain,range,BandedOperator(-1,1, (i,j) -> 2/(dom.b-dom.a)*d(j-1,domain.λ)))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end

function *(E::Evaluation,domain::Ultraspherical)
    if typeof(E.GD) <: NullGrid
        range = GridValues(domain.GD) # Should adapt based on something?
    else
        range = GridValues(E.GD) # Should adapt based on something?
    end
    op = UltrasphericalEvaluation(domain.λ,range)
    ConcreteOperator(domain,range,op)
end

function *(E::Evaluation,C::ConcreteOperator{Dom,Ultraspherical}) where Dom <: Basis
    λ = C.range.λ
    if typeof(E.GD) <: NullGrid
        range = GridValues(C.domain.GD) # Should adapt based on something?
    else
        range = GridValues(E.GD) # Should adapt based on something?
    end
    L = UltrasphericalEvaluation(λ,range)*C.L
    ConcreteOperator(C.domain,range,L)
end

function *(D::Derivative,Dc::ConcreteOperator{Dom,Ultraspherical}) where Dom <: Basis
    if D.order == 1 
        a = Dc.range.GD.D.a
        b = Dc.range.GD.D.b
        λ = Dc.range.λ
        range = Ultraspherical(λ+1,Dc.range.GD)
        domain = Dc.domain
        L = BandedOperator( -1, 1, (i,j) -> 2/(b-a)*d(j-1,λ))*Dc.L
        ConcreteOperator(domain,range,L)
    else
        Derivative(D.order-1)*(Derivative(1)*Dc)
    end
end

function *(Op::CollocatedOperator,domain::Ultraspherical)
    Evaluation(Op.GD)*(Op.Op*domain)
end

function poly(λ::Number,n,z::Vector)
    a, b = Jacobi_ab(λ - 0.5,λ - 0.5)
    poly(a,b,n,z)
end

function UltrasphericalEvaluation(λ,GV::GridValues)
    a, b = Jacobi_ab(λ - 0.5,λ - 0.5)
    OPEvaluationOperator(GV.GD.grid,a,b)
end