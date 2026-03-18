# API Reference

## Module

```@docs
OperatorApproximation
```

## BasisExpansion

```@autodocs
Modules = [OperatorApproximation]
Filter = t -> t isa Type && t <: OperatorApproximation.BasisExpansion || false
```

### Functions

```@docs
BasisExpansion
pad
chop
roots
testconv
moment
setbasis
setN
setgrid
```

## Domains

```@docs
Domain
Interval
Circle
Axis
SemiAxis
ChebyshevMappedInterval
LobattoMappedInterval
UltraMappedInterval
JacobiMappedInterval
PeriodicMappedInterval
UnitCircle
HermiteRealAxis
LaguerreSemiAxis
RationalRealAxis
arclength
isin
```

## Bases

```@docs
Ultraspherical
Jacobi
Fourier
Laurent
Hardy
HermitePoly
HermiteFun
LaguerrePoly
LaguerreFun
GridValues
FixedGridValues
OscRational
MarchenkoPastur
```

## Abstract Operators

```@docs
Derivative
FloquetDerivative
Multiplication
FastMultiplication
Conversion
FastConversion
CauchyTransform
CauchyOperator
BoundaryValue
Residue
Truncation
```

## Operator Algebra

```@docs
diagm
matrix2BlockOperator
```

## Solvers

```@docs
eigen
ContinuousEigen
```

## Riemann–Hilbert

```@docs
RHP
RHSolver
RHSolverVec
rhdomain
rhrange
rhmult
rhrhs
rhplot
rhsplot
```

## Plotting

```@docs
coefplot
domainplot
weightplot
```

## Utilities

```@docs
combine
simp
e
Chop
lancz
RecCoef
clearCauchycache
```
