function isconvertible(b1::DiscreteBasis,b2::Ultraspherical)
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: UltraInterval && b1.GD.λ ≈ b2.λ) || 
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: UltraMappedInterval && b1.GD.λ ≈ b2.λ)
end

function isconvertible(b1::DiscreteBasis,b2::Jacobi)
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: JacobiInterval && b1.GD.α ≈ b2.α && b1.GD.β ≈ b2.β) || 
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: JacobiMappedInterval && b1.GD.α ≈ b2.α && b1.GD.β ≈ b2.β)
end

function isconvertible(b1::DiscreteBasis,b2::Fourier)
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: PeriodicInterval) || 
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: PeriodicMappedInterval)
end

function isconvertible(b1::DiscreteBasis,b2::Rational)
    iscompatible(b1.GD,b2.GD)
end

function isconvertible(b1::DiscreteBasis,b2::Laurent)
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: PeriodicCircle) || 
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: PeriodicMappedCircle)
end

function isconvertible(b1::DiscreteBasis,b2::Hermite)
    iscompatible(b1.GD,b2.GD)
end

function conversion(b1::GridValues,b2::Fourier)
    Op = DiscreteFourierTransform()
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::GridValues,b2::Rational)
    Op = DiscreteFourierTransformII()
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::GridValues,b2::Ultraspherical)
    λ = b2.λ
    a,b = Jacobi_ab(λ - 1/2, λ - 1/2)
    Op = OPEigenTransform(a,b)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::GridValues,b2::HermitePoly)
    a,b = Hermite_ab()
    Op = OPEigenTransform(a,b)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::GridValues,b2::HermiteFun)
    a,b = Hermite_ab()
    Op = OPWeightedEigenTransform(a,b,x -> (2*pi)^(.25)*exp.(x.^2/4))
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::GridValues,b2::Jacobi)
    α = b2.α
    β = b2.β
    a,b = Jacobi_ab(α,β)
    Op = OPEigenTransform(a,b)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::GridValues,b2::Laurent)
    Op = DiscreteFourierTransform()
    ConcreteOperator(b1,b2,Op)
end