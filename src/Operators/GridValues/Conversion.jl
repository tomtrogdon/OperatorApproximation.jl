function isconvertible(b1::DiscreteBasis,b2::Ultraspherical)
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: UltraInterval && b1.GD.λ ≈ b2.λ) || 
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: UltraMappedInterval && b1.GD.λ ≈ b2.λ)
end

function isconvertible(b1::DiscreteBasis,b2::Fourier)
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: PeriodicInterval) || 
    (iscompatible(b1.GD,b2.GD) && typeof(b1.GD) <: PeriodicMappedInterval)
end

function conversion(b1::GridValues,b2::Fourier)
    Op = DiscreteFourierTransform()
    ConcreteLazyOperator(b1,b2,Op)
end

function conversion(b1::GridValues,b2::Ultraspherical)
    λ = b2.λ
    a,b = Jacobi_ab(λ - 1/2, λ - 1/2)
    Op = OPEigenTransform(a,b)
    ConcreteLazyOperator(b1,b2,Op)
end