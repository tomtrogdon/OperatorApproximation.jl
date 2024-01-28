struct Ultraspherical <: Basis
    λ::Number
    GD::GridDomain
end

function (P::BasisExpansion{Ultraspherical})(X::Number) # Clenshaw's algorithm
    n = P.c |> length
    λ = P.basis.λ
    x = P.basis.GD.D.imap(X)
    a,b = Jacobi_ab(λ - 1/2, λ- 1/2)
    (hcat(e(1,n) |> sparse,(Jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
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
        ConcreteLazyOperator(domain,range,SingleBandedOperator(-1,1, (i,j) -> 2/(dom.b-dom.a)*d(j-1,domain.λ)))
    else
        Derivative(D.order-1)*(Derivative(1)*domain)
    end
end

function *(E::Evaluation,domain::Ultraspherical)
    range = GridValues(domain.GD) # inherit the GD
    op = UltrasphericalEvaluation(domain.λ,range)
    ConcreteLazyOperator(domain,range,op)
end

function *(E::Evaluation,C::ConcreteLazyOperator{Dom,Ultraspherical}) where Dom <: Basis
    λ = C.range.λ
    range = GridValues(C.range.GD) # Should adapt based on something?
    L = UltrasphericalEvaluation(λ,range)*C.L
    ConcreteLazyOperator(C.domain,range,L)
end

function *(D::Derivative,Dc::ConcreteLazyOperator{Dom,Ultraspherical}) where Dom <: Basis
    if D.order == 1 
        a = Dc.range.GD.D.a
        b = Dc.range.GD.D.b
        λ = Dc.range.λ
        range = Ultraspherical(λ+1,Dc.range.GD)
        domain = Dc.domain
        L = SingleBandedOperator( -1, 1, (i,j) -> 2/(b-a)*d(j-1,λ))*Dc.L
        ConcreteLazyOperator(domain,range,L)
    else
        Derivative(D.order-1)*(Derivative(1)*Dc)
    end
end

function *(Op::CollocatedOperator,domain::Ultraspherical)
    Evaluation()*(Op.Op*domain)
end

function poly(λ::Number,n,z::Vector)
    a, b = Jacobi_ab(λ - 0.5,λ - 0.5)
    poly(a,b,n,z)
end

function UltrasphericalEvaluation(λ,GV::GridValues)
    a, b = Jacobi_ab(λ - 0.5,λ - 0.5)
    OPEvaluationOperator(GV.GD.grid,a,b)
end

struct UltrasphericalPointFunctional  #TODO: Use parametric type?
    basis::Ultraspherical
    c::Number
end

function *(Lf::LeftBoundaryFunctional,basis::Ultraspherical)
    UltrasphericalPointFunctional(basis,basis.GD.D.a)
end

function *(Lf::RightBoundaryFunctional,basis::Ultraspherical)
    UltrasphericalPointFunctional(basis,basis.GD.D.b)
end

function Matrix(UPF::UltrasphericalPointFunctional,k,n)
    c = UPF.basis.GD.D.imap(UPF.c)
    λ = UPF.basis.λ
    V = poly(λ,n,[c])
    for i = 1:k-1
        V = vcat(V,poly(λ + i,n,[c])*Matrix(Derivative(i)*UPF.basis,n,n))
    end
    V
end

function Matrix(BF::ConcreteBoundaryFunctional{T},n) where T <: Ultraspherical
    Lr = RightBoundaryFunctional()*BF.domain
    Ll = LeftBoundaryFunctional()*BF.domain
    k = size(BF.A)[1]
    BF.A*Matrix(Ll,k,n) + BF.B*Matrix(Lr,k,n)
end
