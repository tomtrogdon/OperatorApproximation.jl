# Empty abstract type to avoid double definitions
abstract type Empty_Abst end

## Duplicate exports: eigen, MappedInterval

## Not Implemented Functions: GridInterval, FiniteGridValues, Transform, setgrid, rhmult, rhrhs, RHSolverVec, rhsplot, combinebasexp

## TO DO: Conversion

## Not Understood by Wietse 😲: HermiteFun, Hardy (what o means), CauchyOperator (what o means), BoundaryValue, RHSolver, RHP, adapt, dilog, rhwellposed, n_partial_horner, p_partial_horner

## Problem: mfft

"""
    Domain

Abstract type of domain structures and abstract types

# Structs
- `Interval`
- `Circle`
- `Axis`
- `SemiAxis`
- `DiscreteDomain`
"""
function Domain(x::Empty_Abst)
end

"""
    Operator

Abstract type of operator abstract types

# Abstract types
- `AbstractOperator`
"""
function Operator(x::Empty_Abst)
end

"""
    BasisEvaluationOperator <: DenseOperator

Abstract type of basis evaluation operators

# Structs
- `OPEvaluationOperator`
- `WeightedOPEvaluationOperator`
- `GenericEvaluationOperator`
- `FourierEvaluationOperator`
- `LaurentEvaluationOperator`
- `PosLaurentEvaluationOperator`
- `NegLaurentEvaluationOperator`
- `RationalEvaluationOperator`
- `OPCauchyEvaluationOperator`
- `PoleResCauchyEvaluationOperator`
"""
function BasisEvaluationOperator(x::Empty_Abst)
end

"""
    MatrixOperator

Abstract type of matrix operators

# Structs
- `DenseOperator`
"""
function MatrixOperator(x::Empty_Abst)
end

"""
    DenseOperator <: MatrixOperator

Abstract type of dense operators

# Structs
- `BasisEvaluationOperator`
- `NaiveTransform`
- `FastTransform`
- `FixedMatrix`
- `DenseTimesBanded`
- `DenseTimesDense`
- `BandedTimesDense`
- `InverseBasicBandedOperator`
- `GridMultiplication`
"""
function DenseOperator(x::Empty_Abst)
end



"""
    NaiveTransform <: DenseOperator

Abstract type of naive transforms

# Structs
- `OPEigenTransform`
- `OPWeightedEigenTransform`
"""
function NaiveTransform(x::Empty_Abst)
end

"""
    FastTransform <: DenseOperator

Abstract type of fast transforms

# Structs
- `DiscreteFourierTransform`
- `DiscreteFourierTransformII`
- `DiscreteCosineTransform`
- `IDiscreteCosineTransform`
"""
function FastTransform(x::Empty_Abst)
end

"""
    GridDomain

Abstract type of grid domain structures and abstract type

# Structs
- `Grid`
- `GridRegion`
- `GridInterval`
- `GridAxis`
- `GridSemiAxis`
- `GridCircle`
"""
function GridDomain(x::Empty_Abst)
end

"""
    Basis

Abstract type of basis structures

# Structs
- `FiniteBasis`
- `AnyBasis`
- `DirectSum`
- `Fourier`
- `DiscreteBasis`
- `Hardy{T <: GridRegion, S <: Domain}`
- `Hermite`
- `Jacobi`
- `Laguerre`
- `Laurent`
- `MarchenkoPastur`
- `OscRational`
- `Ultraspherical`
"""
function Basis(x::Empty_Abst)
end 


"""
    AbstractOperator <: Operator

Abstract type of operators

# Operators    
- `FastConversion` 
- `Multiplication`
- `Conversion`
- `FloquetDerivative`
- `CauchyTransform`
- `BoundaryValue`
- `CauchyOperator`
- `BlockDiagonalAbstractOperator`
- `AbstractZeroOperator`
- `BlockAbstractOperator`
- `CoefConversion`
- `Residue`
- `Truncation`
- `FourierTransform(o::Int64)`
- `Shift`
"""
function AbstractOperator(x::Empty_Abst)
end 

"""
    Derivative(order::Integer)

'order'th order derivative structure.

---

    Derivative()

First order derivative structure.
"""
function Derivative(x::Empty_Abst)
end

"""
    GridInterval

Not implemented
"""
function GridInterval(x::Empty_Abst)
end

"""
    Ultraspherical(λ::Number, GD::GridInterval)   <: Basis   

Ultraspherical (Gegenbauer) type structure with ultrasphericity parameter `λ` and `GridInterval` (`GridDomain`) `GD`.   

# Examples
    julia> λ =  0.0;
    julia> GD = UltraMappedInterval(0.,2.,2.0);
    julia> Ultraspherical(λ,GD);
"""
function Ultraspherical(x::Empty_Abst)
end

"""
    ChebyshevInterval(D::Interval, grid::Function) <: GridInterval

Chebyshev interval structure over `interval` `D` with coefficient number -> grid values function `grid`. 
---

    ChebyshevInterval()

Chebyshev interval structure over [-1,1]
"""
function ChebyshevInterval(x::Empty_Abst)
end

"""
    GridValues(GD::T) {T <: GridDomain} <: DiscreteBasis

GridValues structure
"""
function GridValues(x::Empty_Abst)
end

"""
    FixedGridValues(pts::Vector, GD::GridDomain) <: DiscreteBasis

FixedGridValues structure
"""
function FixedGridValues(x::Empty_Abst)
end

"""
    FiniteGridValues

Not implemented
"""
function FiniteGridValues(x::Empty_Abst)
end

"""
    ConcreteOperator(domain::D, range::R, L::T){D<:Basis,R<:Basis,T<:MatrixOperator}

ConcreteOperator structure of operator `MatrixOperator` `L` with a `Basis` domain `domain` and `Basis` range `range`.
"""
function ConcreteOperator(x::Empty_Abst)
end

"""
     Multiplication(f::Union{Function,BasisExpansion,Vector}) <: AbstractOperator
    
Multiplication structure by `f`.
"""
function Multiplication(x::Empty_Abst)
end

"""
    FastConversion(range::Basis) <: AbstractOperator
    
fast conversion structure to `Basis` `range`

---

    fastconversion(b1::Ultraspherical, b2::GridValues{T}) (T <: Union{UltraMappedInterval,UltraInterval})
    
Constructs `FastConversion` structure that converts `Ultraspherical` structure `b1` to domain of `b2`
"""
function FastConversion(x::Empty_Abst)
end

"""
    ChebyshevMappedInterval(D::Interval, grid::Function) <: GridInterval

Chebyshev mapped interval structure over `interval` `D` with coefficient number -> grid values function `grid`. 

---
    ChebyshevMappedInterval(a,b)

Constructs 'ChebyshevMappedInterval' structure over [`a`, `b`]. 

"""
function ChebyshevMappedInterval(x::Empty_Abst)
end

"""
    MappedInterval(map::Function, imap::Function, a, b) <: Interval

Mapped interval structure that creates a map `map` from [-1,1] to [`a`,`b`] with inverse map `imap`.

---
    MappedInterval(a,b)

Constructs `MappedInterval` structure for [`a`, `b`].
"""
function MappedInterval(x::Empty_Abst)
end

"""
    BasisExpansion(basis::T, c::Vector) {T<:Basis} 

Basis expansion structure of `Basis` `basis` with coefficient vector `c`.

---

    BasisExpansion(f::Function,basis::Basis,N::Integer)

Constructs `BasisExpansion` structure for function `f` over `Basis` `basis` using `N` coefficients `c`

# Example

    julia> gd = UltraMappedInterval(0., 2., 0.);
    julia> sp = Ultraspherical(2.0, gd);    
    julia> f =  BasisExpansion(x -> sin(x), sp, 24);
---

    BasisExpansion(f::Function, basis::Basis)

Constructs `BasisExpansion` structure for function `f` over `Basis` `basis`.
The amount of coefficients are adaptively decided based on accuracy.

# Example

    julia> gd = UltraMappedInterval(0.,2.,0.);
    julia> sp = Ultraspherical(2.0, gd);    
    julia> f =  BasisExpansion(x -> sin(x), sp);
---

    BasisExpansion(f::BasisExpansion, sp::Basis)

Converts `BasisExpansion` `f` to `Basis` `sp`

---

    BasisExpansion(f::BasisExpansion, sp::Basis, N::Integer)

Converts `BasisExpansion` `f` to `Basis` `sp` with `N` total coefficients through padding

---

    BasisExpansion(f::Function, basis::GridValues, N::Integer)

Constructs `BasisExpansion` structure using basis `basis` and `N` coefficients c, the value of `f` at the grid values. 
"""
function BasisExpansion(x::Empty_Abst)
end

"""
    Conversion(range::Basis) <: AbstractOperator
    
Conversion structure 

# TO COMPLETE

"""
function Conversion(x::Empty_Abst)
end

"""
Laurent

    conversion(b1::Laurent,b2::GridValues)
    
    conversion(b1::Laurent,b2::FixedGridValues)

    conversion(b1::Laurent,b2::Laurent)

    conversion(b1::GridValues,b2::Laurent)

---
OscRational

    conversion(b1::OscRational,b2::OscRational)
     
    conversion(b1::OscRational,b2::GridValues)

    conversion(b1::OscRational,b2::FixedGridValues)

    conversion(b1::GridValues,b2::OscRational)

---
Fourier

    conversion(b1::Fourier,b2::GridValues)

    conversion(b1::Fourier,b2::FixedGridValues)

    conversion(b1::Fourier,b2::Fourier)

---
Ultraspherical

    conversion(b1::GridValues,b2::Ultraspherical)

    conversion(b1::Ultraspherical,b2::GridValues{T}) where T

    conversion(b1::Ultraspherical,b2::FixedGridValues)

    conversion(b1::Ultraspherical,b2::Ultraspherical)

---
MarchenkoPastur

    conversion(b1::GridValues,b2::MarchenkoPastur)

---
HermitePoly

    conversion(b1::GridValues,b2::HermitePoly)
    
    conversion(b1::HermitePoly,b2::GridValues)

    conversion(b1::HermitePoly,b2::FixedGridValues)

    conversion(b1::HermitePoly,b2::HermitePoly)


---
HermiteFun

    conversion(b1::GridValues,b2::HermiteFun)

    conversion(b1::HermiteFun,b2::HermiteFun)

    conversion(b1::HermiteFun,b2::GridValues)

    conversion(b1::HermiteFun,b2::FixedGridValues)

---
LaguerrePoly

    conversion(b1::GridValues,b2::LaguerrePoly)

--- 
LaguerreFun

    conversion(b1::GridValues,b2::LaguerreFun)

---
Jacobi

    conversion(b1::GridValues,b2::Jacobi)

    conversion(b1::Jacobi,b2::GridValues)

    conversion(b1::Jacobi,b2::FixedGridValues)

---
Hardy

    conversion(b1::Hardy{Exterior{T},S},b2::GridValues) where {T <: Union{JacobiMappedInterval,JacobiInterval}, S <: Interval}

    conversion(b1::Hardy{Exterior{T},S},b2::GridValues) where {T <: Union{MarchenkoPasturMappedInterval,MarchenkoPasturInterval}, S <: Interval}

    conversion(b1::Hardy{Exterior{T},S},b2::FixedGridValues) where {T <: Union{JacobiMappedInterval,JacobiInterval}, S <: Interval}

    conversion(b1::Hardy{Exterior{T},S},b2::FixedGridValues) where {T <: Union{MarchenkoPasturMappedInterval,MarchenkoPasturInterval}, S <: Interval}

    conversion(b1::Hardy{T,S},b2::GridValues) where {T <: Exterior, S <: DiscreteDomain}

    conversion(b1::Hardy{T,S},b2::FixedGridValues) where {T <: Exterior, S <: DiscreteDomain}

Conversion functions

"""
function conversion(x::Empty_Abst)
end

"""
    UnitInterval()

Construct unit interval structure `UnitInterval`.
"""
function UnitInterval(x::Empty_Abst)
end

"""
    Transform

Not implemented
"""
function Transform(x::Empty_Abst)
end

"""
    setbasis(b::Basis)

Sets global basis `b`
"""
function setbasis(x::Empty_Abst)
end

"""
    setgrid

Not implemented
"""
function setgrid(x::Empty_Abst)
end

"""
    setN(N)

Sets global coefficient number `N`. Either integer or "adaptive"
"""
function setN(x::Empty_Abst)
end

"""
    UltraInterval(D::Interval, λ::Number, grid::Function) <: GridInterval

Ultraspherical interval structure with `Interval` `D`, ultrasphericity `λ` and coefficient number -> grid values function `grid`.
Typically over interval [-1,1].

---

    UltraInterval(λ)

Constructs `UltraInterval` with ultrasphericity `λ` over [-1,1] with corresponding Gauss quadrature grid.
"""
function UltraInterval(x::Empty_Abst)
end

"""
    JacobiInterval(D::Interval, α::Number, β::Number, grid::Function) <: GridInterval

Jacobi polynomial interval structure with `Interval` `D`, weight function (1-x)`ᵅ`(1+x)`ᵝ` and coefficient number -> grid values function `grid`.
Typically over interval [-1,1].

---

    JacobiInterval(α, β)

Constructs `JacobiInterval` with weight function (1-x)`ᵅ`(1+x)`ᵝ` over [-1,1] with corresponding Gauss quadrature grid.

"""
function JacobiInterval(x::Empty_Abst)
end

"""
    UltraMappedInterval(D::Interval, λ::Float64, grid::Function) <: GridInterval
    
Mapped Ultraspherical interval structure with `Interval` `D`, ultrasphericity `λ` and coefficient number -> grid values function `grid`.

---

    UltraMappedInterval(a, b, λ)

Constructs `UltraMappedInterval` over interval [`a`, `b`] with ultrasphericity `λ` with corresponding Gauss quadrature grid.
    
"""
function UltraMappedInterval(x::Empty_Abst)
end

"""
    JacobiMappedInterval(D::Interval, α::Number, β::Number, grid::Function) <: GridInterval

Mapped Jacobi polynomial interval structure with `Interval` `D`, weight function (1-x)`ᵅ`(1+x)`ᵝ` and coefficient number -> grid values function `grid`.

---

    JacobiMappedInterval(a, b, α, β)

Constructs `JacobiMappedInterval` over interval [`a`, `b`] with weight function (1-x)`ᵅ`(1+x)`ᵝ` with corresponding Gauss quadrature grid.

"""
function JacobiMappedInterval(x::Empty_Abst)
end

"""
    PeriodicInterval(D::Interval, grid::Function) <: GridInterval

Periodic interval structure with `Interval` `D` and coefficient number -> grid values function `grid`.

---

    PeriodicInterval()

Construct `PeriodicInterval` over interval [-1, 1] with uniform grid mapping.

"""
function PeriodicInterval(x::Empty_Abst)
end

"""
    PeriodicMappedInterval(D::Interval, grid::Function) <: GridInterval

Mapped periodic interval structure with `Interval` `D` and coefficient number -> grid values function `grid`.

---

    PeriodicMappedInterval(a,b)

Construct `PeriodicInterval` over interval [`a`, `b`] with uniform grid mapping.

"""
function PeriodicMappedInterval(x::Empty_Abst)
end

"""
    Fourier(GD::GridInterval) <: Basis

Fourier structure
"""
function Fourier(x::Empty_Abst)
end

"""
    ⊞(A1::{AbstractOperator,BlockAbstractOperator}, A2::{AbstractOperator,BlockAbstractOperator}, ...)

Horizontaly concatenates `AbstractOperator`s and/or `BlockAbstractOperator`s as `BlockAbstractOperator`.

---

    ⊞(A1::{MatrixOperator,BlockMatrixOperator},A2::{MatrixOperator,BlockMatrixOperator})

Horizontaly concatenates `MatrixOperator`s and/or `BlockMatrixOperator`s as `BlockMatrixOperator`.

---

    ⊞(A1::ConcreteOperator,A2::ConcreteMatrixOperator)

Horizontaly concatenates `ConcreteOperator`s as a `ConcreteOperator`.    
"""
function ⊞(x::Empty_Abst)
end

"""
    DirectSum(bases::Vector{T}) {T <: Basis} <: Basis

Direct sum structure of bases vector `Basis`

---

    DirectSum(b::Union{Vector{S},S}) where S <: Basis

Constructs `DirectSum` of bases and/or bases vectors.
"""
function DirectSum(x::Empty_Abst)
end

"""
    ⊘(A1::{AbstractOperator,BlockAbstractOperator}, A2::{AbstractOperator,BlockAbstractOperator}, ...)

Vertically concatenates `AbstractOperator`s and/or `BlockAbstractOperator`s as `BlockAbstractOperator`.

---

    ⊘(A1::{MatrixOperator,BlockMatrixOperator},A2::{MatrixOperator,BlockMatrixOperator})

Vertically concatenates `MatrixOperator`s and/or `BlockMatrixOperator`s as `BlockMatrixOperator`.

---

    ⊘(A1::ConcreteOperator,A2::ConcreteMatrixOperator)

Vertically concatenates `ConcreteOperator`s as a `ConcreteOperator`.    
"""
function ⊘(x::Empty_Abst)
end

"""
    ⊕(b1::{Basis,DirectSum},b2::{Basis,DirectSum})

Combines `b1` and `b2` in a `DirectSum` structure.

---

    ⊕(f::BasisExpansion{T},g::BasisExpansion{S})

Combines `f` and `g` in a `BasisExpansion` structure with `Basis` `DirectSum`.

---

    ⊕(A1::{AbstractOperator, BlockDiagonalAbstractOperator},A2::{AbstractOperator, BlockDiagonalAbstractOperator})

Combines `A1` and `A2` in a `BlockDiagonalAbstractOperator` structure.

"""
function ⊕(x::Empty_Abst)
end

"""
    FloquetDerivative(order::Integer, μ::Float64) <: AbstractOperator

`order`th order Floquet derivative structure with Floquet multipliers `μ`

---

FloquetDerivative(μ::Float64)

First order Floquet derivative constructor with Floquet multipliers `μ`
"""
function FloquetDerivative(x::Empty_Abst)
end

"""
    ploteval(E::ContinuousEigen)

Scatters eigenvalues of `E`
"""
function ploteval(x::Empty_Abst)
end

"""
    UnitCircle(map::Function, imap::Function, cen::ComplexF64, rad::Float64) <: Circle  
    
Unit circle structure with center `cen`, radius `rad`, and mapping `map` and inverse mapping `imap`
    
--- 

    UnitCircle()

Constructs the unit circle with center 0 and radius 1 in `UnitCircle` structure.
    
"""
function UnitCircle(x::Empty_Abst)
end

"""
    MappedCircle(map::Function, imap::Function, cen::ComplexF64, rad::Float64) <: Circle  
    
Unit circle structure with center `cen`, radius `rad`, and mapping `map` and inverse mapping `imap`
    
--- 

    UnitCircle(cen,rad)

Constructs the unit circle with center `cen` and radius `rad` in `UnitCircle` structure.
    
"""
function MappedCircle(x::Empty_Abst)
end

"""
    PeriodicCircle(D::Circle, grid::Function) <: GridCircle
    
Periodic circle structure with 'Circle' domain 'D' and coefficient number -> grid values function `grid`

---

    PeriodicCircle()

Unit circle with uniform gridpoints in `PeriodicCircle` structure.
"""
function PeriodicCircle(x::Empty_Abst)
end

"""
    PeriodicMappedCircle(D::Circle, grid::Function) <: GridCircle
    
Periodic circle structure with 'Circle' domain 'D' and coefficient number -> grid values function `grid`

---

    PeriodicMappedCircle(cen,rad)

Circle with center `cen` and radius `rad` with uniform gridpoints in `PeriodicMappedCircle` structure.
"""
function PeriodicMappedCircle(x::Empty_Abst)
end

"""
    Laurent(GD::GridCircle) <: Basis

Laurent polynomials structure on `GridCircle` `GD` (`GridDomain`).
    
"""
function Laurent(x::Empty_Abst)
end

"""
    struct Hardy{T <: GridRegion, S <: Domain}(GD::T) <: Basis

Hardy structure

---

    Hardy(GD) = Hardy{typeof(GD),typeof(GD.GD.D)}(GD)

"""
function Hardy(x::Empty_Abst)
end

"""
    Jacobi(α::Number, β::Number, GD::GridInterval) <: Basis
    
Jacobi polynomial structure using weight function (1-x)`ᵅ`(1+x)`ᵝ` over grid interval `GD`.
    
"""
function Jacobi(x::Empty_Abst)
end

"""
    CauchyTransform <: AbstractOperator

Cauchy transform structure under `AbstractOperator`.

"""
function CauchyTransform(x::Empty_Abst)
end

"""
    Exterior{T <: GridDomain}(D::Domain, grid::Union{Function,Vector}, GD::T) <: GridRegion
    
Exterior structure

---

    Exterior{T}(GD::T) where T <: GridDomain
    Exterior(GD)

---

Constructs `Exterior` structure from `GD`.
"""
function Exterior(x::Empty_Abst)
end

"""
    Interior{T <: GridDomain}(D::Domain, grid::Union{Function,Vector}, GD::T) <: GridRegion
    
Interior structure

---

    Interior{T}(GD::T) where T <: GridDomain
    Interior(GD)

---

Constructs `Interior` structure from `GD`.
"""
function Interior(x::Empty_Abst)
end

"""
    CauchyOperator(o::Int64) <: AbstractOperator
    
Cauchy operator structure with parameter `o`

"""
function CauchyOperator(x::Empty_Abst)
end

"""
    ArgNum(z::ComplexF64, ρ::Float64, θ::Float64) <: Number
    
Complex number in polar structure: `z`=`ρ`e`ᶿ` 

---

    ArgNum(z::Number,ρ::Float64,θ::Float64)

Constructs `ArgNum` of `z`. Converts `z` to type `ComplexF64` and `θ` ∈ [0,2π)
    
"""
function ArgNum(x::Empty_Abst)
end

"""
    LobattoMappedInterval(D::Interval, grid::Function) <: GridInterval

Lobatto mapped interval structure over `interval` `D` with coefficient number -> grid values function `grid`. 

---

    LobattoMappedInterval(a,b)

Constructs `LobattoMappedInterval` over interval [`a`,`b`] with Gauss-Lobatto grid.
"""
function LobattoMappedInterval(x::Empty_Abst)
end

"""
    LobattoInterval(D::Interval, grid::Function) <: GridInterval

Lobatto interval structure over `interval` `D` with coefficient number -> grid values function `grid`. 

---

    LobattoInterval(a,b)

Constructs `LobattoInterval` over interval [-1,1] with Gauss-Lobatto grid.
"""
function LobattoInterval(x::Empty_Abst)
end

"""
    BoundaryValue(o::Int64, range::Basis) <: AbstractOperator
    
Boundary value structure

---

    BoundaryValue(o::Int64,ran::DirectSum)

Constructs `BlockDiagonalAbstractOperator` using `BoundaryValue` as function and the bases of `ran`.
"""
function BoundaryValue(x::Empty_Abst)
end

"""
    BlockDiagonalAbstractOperator{T}(Ops::Vector{T}) <: AbstractOperator where T <: AbstractOperator
    
Vector of block diagonal operators structure. 
"""
function BlockDiagonalAbstractOperator(x::Empty_Abst)
end

"""
    struct AbstractZeroOperator <: AbstractOperator 

Zero operator structure.

---

    AbstractZeroOperator(n::Int64,m::Int64) 

Constructs `m`×`n` `BlockDiagonalAbstractOperator` of `AbstractZeroOperator`.
"""
function AbstractZeroOperator(x::Empty_Abst)
end

"""
    ZeroOperator{T<:CoefficientDomain, S<: CoefficientDomain}(nm::Integer, np::Integer, A::Function) <: SingleBandedOperator

Zero operator structure 

---

    ZeroOperator() = ZeroOperator{𝕏,𝕏}(0,0,x -> 0)
    
Constructs `ZeroOperator` on {𝕏,𝕏}.
"""
function ZeroOperator(x::Empty_Abst)
end

"""
    DirectedLobattoMappedInterval(D::Interval, grid::Function, dgrid::Function) <: DirectedGridInterval
    
Directed Lobatto mapped interval over `Interval` `D` with coefficient number -> grid values function `grid` and directed version `dgrid`.
    
---

    DirectedLobattoMappedInterval(a,b)

Constructs `DirectedLobattoMappedInterval` over interval [`a`, `b`] with directed Gauss-Lobatto grid.
"""
function DirectedLobattoMappedInterval(x::Empty_Abst)
end

"""
    DirectedLLobattoMappedInterval(D::Interval, grid::Function, dgrid::Function) <: DirectedGridInterval
    
Left directed Lobatto mapped interval over `Interval` `D` with coefficient number -> grid values function `grid` and directed version `dgrid`.
    
---

    DirectedLLobattoMappedInterval(a,b)

Constructs `DirectedLLobattoMappedInterval` over interval [`a`, `b`] with directed Gauss-Lobatto grid without `b`.
"""
function DirectedLLobattoMappedInterval(x::Empty_Abst)
end

"""
    DirectedRLobattoMappedInterval(D::Interval, grid::Function, dgrid::Function) <: DirectedGridInterval
    
Right directed Lobatto mapped interval over `Interval` `D` with coefficient number -> grid values function `grid` and directed version `dgrid`.
    
---

    DirectedRLobattoMappedInterval(a,b)

Constructs `DirectedRLobattoMappedInterval` over interval [`a`, `b`] with directed Gauss-Lobatto grid without `a`.
"""
function DirectedRLobattoMappedInterval(x::Empty_Abst)
end

"""
    Legendre(a,b)

Constructs Legendre polynomial structure as `Jacobi` structure with α, β = 0 over interval [`a`, `b`].
"""
function Legendre(x::Empty_Abst)
end

"""
    rhrange(D::Basis)

Constructs `GridValues` structure of `DirectedLobattoMappedInterval` with the interval of `D`.

---

    rhrange(D::DirectSum)

Constructs `Directsum` structure of `GridValues` of `DirectedLobattoMappedInterval` with the interval of `Basis` of `D`.

"""
function rhrange(x::Empty_Abst)
end

"""
    rhdomain(A::Matrix)

Constructs `DirectSum` of `Legendre` type structures with endpoints the collumn elements of `A`∈ Rⁿ×²   
"""
function rhdomain(x::Empty_Abst)
end

"""
    BlockAbstractOperator{T}(Ops::Matrix{T}) <: AbstractOperator where T <: AbstractOperator
    
Block abstract operator structure

---

    BlockAbstractOperator(Op::AbstractOperator,n,m)

Constructs `BlockAbstractOperator` of `n`×`m` `Matrix` filled with `Op`.

"""
function BlockAbstractOperator(x::Empty_Abst)
end

"""
    rhmult

Not Implemented
"""
function rhmult(x::Empty_Abst)
end

"""
    rhrhs

Not Implemented
"""
function rhrhs(x::Empty_Abst)
end

"""
    matrix2BlockOperator(M::Matrix{T}) where T <: Operator

Converts `Matrix` if `Operator`s to `BlockOperator` version.
"""
function matrix2BlockOperator(x::Empty_Abst)
end

"""
    RHSolver(S::ConcreteOperator, jumps, res)

Riemann Hilbert solver structure with solution operators `S`, jumps `Jumps` and residues `res`.

---

    RHSolver(S,jumps)

Constructs `RHSolver` with no residues.

---

    RHSolver(rhp::RHP{T}) where T <: {Matrix, Vector}

Constructs `RHSolver` from the `Matrix` or `Vector` Riemann Hilbert problem `rhp`.

--- 

# Operator

    (R::RHSolver)(c,n)

?

"""
function RHSolver(x::Empty_Abst)
end

"""
    domainplot(T::{Interval, Circle, DiscreteDomain, BasisExpansion, Basis}; kwargs...) 

Plots the domain `T`.
"""
function domainplot(x::Empty_Abst)
end

"""
    mvf2mof(f,n,m)

Lazy conversion of function `f` with matrix output into `n`×`m` matrix of functions outputs `fᵢⱼ`.
"""
function mvf2mof(x::Empty_Abst)
end

"""
    coefplot(f::BasisExpansion{T}; kwargs...)

Plots coefficients of `f`.
"""
function coefplot(x::Empty_Abst)
end

"""
    RHSolverVec

Not implemented
"""
function RHSolverVec(x::Empty_Abst)
end

"""
    arclength(I::{Interval,Circle,GridDomain})

Computes the arc length of `I`.
"""
function arclength(x::Empty_Abst)
end

"""
    RHP{T <: Union{Matrix,Vector}}(Γ::T, J::Vector, P::Vector, R::Vector)
    
Riemann Hilbert Problem structure characterized by a contour `Γ`` (`Matrix` or `Vector`), and jump condition `J`, prescribed values `P`, and residues `R`.

---

    RHP(Γ::T, J::Vector)

Constructs `RHP` with no prescripted valued or residues
"""
function RHP(x::Empty_Abst)
end

"""
    adapt(rhp::RHP,j,ϵ::Float64)

???
"""
function adapt(x::Empty_Abst)
end

"""
    mofeval(f,z)

Matrix of functions evaluation.
"""
function mofeval(x::Empty_Abst)
end

"""
    mult2x2(A,B)

2x2 block "matrix" multiplication.
"""
function mult2x2(x::Empty_Abst)
end

"""
    
    dilog(z)

Computes dilogarithm of `z`
"""
function dilog(x::Empty_Abst)
end

"""
    rhwellposed(rhp::RHP)

???
"""
function rhwellposed(x::Empty_Abst)
end

"""
    rhsplot

Not implemented.
"""
function rhsplot(x::Empty_Abst)
end

"""
    rhplot(rhp::RHP;kwargs...)

Plots the Riemann Hilber Problem.
"""
function rhplot(x::Empty_Abst)
end

"""
    clearCauchycache()

Clears Cauchy cache.
"""
function clearCauchycache(x::Empty_Abst)
end

"""
    HermitePoly(GD::GridAxis) <: Hermite

Hermite polynomial structure over axis `GD`.
"""
function HermitePoly(x::Empty_Abst)
end

"""
    HermiteFun(GD::GridAxis) <: Hermite

Hermite function structure over axis `GD`.
"""
function HermiteFun(x::Empty_Abst)
end

"""
    Axis <: Domain

Abstract type of axis structures
- `RealAxis`
- `MappedAxis`
"""
function Axis(x::Empty_Abst)
end

"""
    GridAxis <: GridDomain

Abstract type of axis structures
- `RationalRealAxis`
- `HermiteRealAxis`
"""
function GridAxis(x::Empty_Abst)
end

"""
    RealAxis(map::Function, imap::Function, cen::Union{ComplexF64,Float64}, θ::Float64) <: Axis

Real axis structure centered at `cen` at angle `θ` with mapping `map` and inverse mapping `imap`.

---

    RealAxis()

Creates the real axis center at 0 with no angle as `RealAxis` structure.
"""
function RealAxis(x::Empty_Abst)
end

"""
    HermiteRealAxis(D::Axis, grid::Function) <: GridAxis
    
Hermite real axis structure with `Interval` `D` and coefficient number -> grid values function `grid`.

---

    HermiteRealAxis()

Constructs real axis and grid function associated with the Hermite orthogonal polynomials as `HermiteRealAxis` structure.
"""
function HermiteRealAxis(x::Empty_Abst)
end

"""
    CoefConversion(range::Basis) <: AbstractOperator

Coefficient conversion operator. Applied using *.
"""
function CoefConversion(x::Empty_Abst)
end

"""
    Erf(GD::GridAxis) <: FiniteBasis

Error function structure
"""
function Erf(x::Empty_Abst)
end

"""
    lancz(N,x,w)

Lanczos algorithm. --- https://www.cs.purdue.edu/archives/2002/wxg/codes/lanczos.m

Generates the first n recurrence coefficients ab of the 
corresponding discrete orthogonal polynomials and puts it in `SymTridiagonal`

The script is adapted from the routine RKPW in
W.B. Gragg and W.J. Harrod, ``The numerically stable 
reconstruction of Jacobi matrices from spectral data'', 
Numer. Math. 44 (1984), 317-335.
"""
function lancz(x::Empty_Abst)
end

"""
    RecCoef(const ΛW::Function, const NN::Function, J::SymTridiagonal)

Recursion coefficient mutable structure.

---

    RecCoef(ΛW,NN)

Construct `RecCoef` using `ΛW` and `NN` with `J` a 2x2 identity matrix.
"""
function RecCoef(x::Empty_Abst)
end

"""
    DiscreteDomain(map::Function, imap::Function, pts::Vector{ComplexF64}) <: Domain
    
Discrete domain structure storing given `Complex` points `pts` with mapping `map` and inverse mapping `imap`.

    DiscreteDomain(pts::Vector)

Constructs `DiscreteDomain` using points `pts` and the identity map.
    
"""
function DiscreteDomain(x::Empty_Abst)
end

"""
    Grid(D::DiscreteDomain, grid::Vector) <: GridDomain
    
Grid structure, storing the discrete domain and specific points

---

    Grid(D::DiscreteDomain) = Grid(D,D.pts)
    Grid(V::Vector) = Grid(DiscreteDomain(V),V)

"""
function Grid(x::Empty_Abst)
end

"""
    Residue(range::Basis) <: AbstractOperator
    
Residue structure.
"""
function Residue(x::Empty_Abst)
end

"""
    function moment(f::BasisExpansion{DirectSum},k::Int64)

Computes the sum `k`th moments of the compotents of `f`.

---

    moment(f::BasisExpansion{T},k::Int64)

Computes the `k`th moment of function `f`.
"""
function moment(x::Empty_Abst)
end

"""
    Truncation(Op::AbstractOperator, k::Int64) <: AbstractOperator

Truncation operator sturcture. Applied using *. Truncates operator at `k` entries
"""
function Truncation(x::Empty_Abst)
end

"""
    MarchenkoPasturInterval(D::Interval, d::Number, grid::Function) <: GridInterval
    
Marchenko Pastur interval structure over `Interval` `D` with coefficient number -> grid values function `grid` and scaled support width `d`=4√γ

---

    MarchenkoPasturInterval(d)

Constructs `MarchenkoPasturInterval` over interval [(1 - √`d`)², (1 + √`d`)²] with appropriate Marchenko-Pastur grid function.
"""
function MarchenkoPasturInterval(x::Empty_Abst)
end

"""
    MarchenkoPastur(d::Number, GD::GridInterval) <: Basis

Marchenko Pastur type structure with scaled support width `d` and `GridInterval` (`GridDomain`) `GD`.   

"""
function MarchenkoPastur(x::Empty_Abst)
end

"""
    MarchenkoPasturMappedInterval(D::Interval, d::Number, grid::Function) <: GridInterval
    
Mapped Marchenko Pastur interval structure over `Interval` `D` with coefficient number -> grid values function `grid` and scaled support width `d`=4√γ

---

    MarchenkoPasturMappedInterval(d)

Constructs `MarchenkoPasturMappedInterval` over interval [`a`, `b`] with scaled support width `d` and appropriate Marchenko-Pastur grid function.
"""
function MarchenkoPasturMappedInterval(x::Empty_Abst)
end

"""
    weightplot(f::BasisExpansion; dx = .01, kwargs...)

Plots weighted integrad.
"""
function weightplot(x::Empty_Abst)
end

"""
    RationalRealAxis(D::Axis, grid::Function) <: GridAxis

Rational real axis structure with `Interval` `D` and coefficient number -> grid values function `grid`.

---

    RationalRealAxis()

Constructs real axis and grid function associated with the Rational orthogonal polynomials as `RationalRealAxis` structure.
"""
function RationalRealAxis(x::Empty_Abst)
end

"""
    OscRational(GD::GridAxis, α::Number) <: Basis
    
Oscilatory rational function type structure with power `α` and `GridInterval` (`GridDomain`) `GD`.   
"""
function OscRational(x::Empty_Abst)
end

"""
    ⊙(g::Function,f::BasisExpansion{T}) where T <: {GridValues, Ultraspherical}

Function composition operator: g∘f
"""
function ⊙(x::Empty_Abst)
end

# """
    
# """
# function mfft(x::Empty_Abst)
# end

"""
    sumdot(f::BasisExpansion,g::BasisExpansion)

Innerproduct of `f` and `g`.

--- 

    sumdot(f::BasisExpansion{T},g::BasisExpansion{T}) where T <: DirectSum

summed Innerproduct of `f` and `g`. All permutations of T being a `DirectSum` or not is available.

"""
function sumdot(x::Empty_Abst)
end

"""
    combine(f::BasisExpansion{T}) where T <: DirectSum

Adds all compatible bases in `f` into a `BasisExpansion{DirectSum}`
---

    combine(f::Vector{T}) where T <: BasisExpansion

Combines element wise.
"""
function combine(x::Empty_Abst)
end

"""
    combinebasexp

Not implemented
"""
function combinebasexp(x::Empty_Abst)
end

"""
    simp(f::Vector{T}) where T <: BasisExpansion

Simplifies array of `BasisExpansion`. `Chop`s, `combine`s and then `chop`s.

---

Someone who does way too much for a person they like.

"""
function simp(x::Empty_Abst)
end

"""
    RationalMappedAxis(D::Axis, grid::Function) <: GridAxis
    
Rational mapped axis structure with `Interval` `D` and coefficient number -> grid values function `grid`.

---

    RationalMappedAxis()

Constructs mapped axis and grid function associated with the Rational orthogonal polynomials as `RationalMappedAxis` structure.    
"""
function RationalMappedAxis(x::Empty_Abst)
end

"""
    MappedSemiAxis(map::Function, imap::Function, cen::Number, θ::Float64) <: SemiAxis

Mapped half axis structure with center `cen` and rotation angle `θ` of axis with mapping `map` and inverse mapping `imap`.

    MappedSemiAxis(σ, cen, θ)

Constructs `MappedSemiAxis` with center `cen`, rotation angle `θ` and at mapping rate `σ`.
"""
function MappedSemiAxis(x::Empty_Abst)
end

"""
    LaguerreSemiAxis(D::SemiAxis, grid::Function, α::Float64) <: GridSemiAxis

Laguerre half axis structure with `SemiAxis` `D`, coefficient number -> grid values function `grid` and order `α`.

    LaguerreSemiAxis(D, α)

Constructs Laguerre half axis on `SemiAxis` `D` with order `α` and grid function associated with the Laguerre orthogonal polynomials as `LaguerreSemiAxis` structure.    
"""
function LaguerreSemiAxis(x::Empty_Abst)
end

"""
    LaguerrePoly(GD::GridSemiAxis, α::Float64) <: Laguerre
    
Laguerre polynomial type structure with order `α` and `GridSemiAxis` (`SemiAxis`) `GD`.   

"""
function LaguerrePoly(x::Empty_Abst)
end

"""
    LaguerreFun(GD::GridSemiAxis, α::Float64) <: Laguerre

Laguerre type structure for decomposed functions with order `α` and `GridSemiAxis` (`SemiAxis`) `GD`.   

"""
function LaguerreFun(x::Empty_Abst)
end

"""
    roots(P::BasisExpansion{Ultraspherical}; tol = 1e-13)

Finds roots of `P` with tolerance `tol`.

#Example

    julia> λ =  0.0;
    julia> gd = UltraMappedInterval(0.,2.,2.0);
    julia> sp = Ultraspherical(λ, gd);
    julia> f = BasisExpansion(x -> cos(π*x), sp);
    julia> roots(f, tol = 1e-15)
    2-element Vector{Float64}:
     0.49999999999999944
     1.4999999999999998
"""
function roots(x::Empty_Abst)
end

"""
    FourierTransform(o::Int64) <: AbstractOperator

Fourier transform structure with parameter o.
"""
function FourierTransform(x::Empty_Abst)
end






"""
    p_partial_horner_vec(n,x)

???
"""
function p_partial_horner_vec(x::Empty_Abst)
end

"""
    n_partial_horner_vec(n,x)

Used to evaluate Hardy Basis expansion
"""
function n_partial_horner_vec(x::Empty_Abst)
end

"""
    Shift(order::Integer) <: AbstractOperator

Shift structure with order `order`
"""
function Shift(x::Empty_Abst)
end

"""
    isconvertible(b1::Laurent,b2::DiscreteBasis)

    isconvertible(b1::Laurent,b2::Laurent)
    
    isconvertible(b1::OscRational,b2::DiscreteBasis)

    isconvertible(b1::OscRational,b2::OscRational)

    isconvertible(b1::Fourier,b2::DiscreteBasis)

    isconvertible(b1::Fourier,b2::Fourier)

    isconvertible(b1::DiscreteBasis,b2::Ultraspherical)

    isconvertible(b1::DiscreteBasis,b2::Jacobi)

    isconvertible(b1::DiscreteBasis,b2::MarchenkoPastur)

    isconvertible(b1::DiscreteBasis,b2::Fourier)

    isconvertible(b1::DiscreteBasis,b2::OscRational)

    isconvertible(b1::DiscreteBasis,b2::Laurent)

    isconvertible(b1::DiscreteBasis,b2::Hermite)

    isconvertible(b1::DiscreteBasis,b2::Laguerre)

    isconvertible(b1::Hardy,b2::Basis)

    isconvertible(b1::Hermite,b2::DiscreteBasis)

    isconvertible(b1::Hermite,b2::Hermite)

    isconvertible(b1::Jacobi,b2::DiscreteBasis)

    isconvertible(b1::Ultraspherical,b2::DiscreteBasis)

    isconvertible(b1::Ultraspherical,b2::Ultraspherical)

Checks if it's convertible through iscompatible.
"""
function isconvertible(x::Empty_Abst)
end


"""
    toeplitz(c::Vector,N::Integer)

Sparse diagonal Toeplitz matrix with vector `c` and `N` ?
"""
function toeplitz(x::Empty_Abst)
end

"""
    toeplitz_function(c::Vector)

Returns toeplitz operator F(k,J) acting on vector `c`
"""
function toeplitz_function(x::Empty_Abst)
end

"""
    poleres_cauchy(ps,z::Number)

pole residue of `ps` at number `z` (first order poles)

---
    poleres_cauchy(ps,z::Vector)

pole residue of `ps` at vector `z` (first order poles)
"""
function poleres_cauchy(x::Empty_Abst)
end


"""
    herm_d(i,j)

Hermite derivative operator
"""
function herm_d(x::Empty_Abst)
end

"""
    hasfastconversion(b1::Ultraspherical, b2::DiscreteBasis)

Checks if there is fast conversion between the two basis through mod(b1.λ - b2.GD.λ,1) < 1e-15

---

    hasfastconversion(b1::Ultraspherical,b2::GridValues{T}) where T <: Union{UltraMappedInterval,UltraInterval}

Checks if there is fast conversion for an ultrspherical basis (GridValues is from Chebyshev)

"""
function hasfastconversion(x::Empty_Abst)
end

"""
    _s(j,λ)

That one ultraspherical number
"""
function _s(x::Empty_Abst)
end

"""
    _t(j,λ)

That other ultraspherical number
"""
function _t(x::Empty_Abst)
end


"""
    _conv_ultra(λ)

Gives operator that converts between ultraspherical polynomials
"""
function _conv_ultra(x::Empty_Abst)
end