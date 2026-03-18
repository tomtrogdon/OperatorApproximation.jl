# OperatorApproximation.jl

[![CI](https://github.com/tomtrogdon/OperatorApproximation.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/tomtrogdon/OperatorApproximation.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/tomtrogdon/OperatorApproximation.jl/graph/badge.svg?token=QXMLQ083L3)](https://codecov.io/gh/tomtrogdon/OperatorApproximation.jl)
[![DOI](https://zenodo.org/badge/742642635.svg)](https://zenodo.org/doi/10.5281/zenodo.10957640)

**OperatorApproximation.jl** is a Julia framework for approximating functions and operators, and for solving operator equations using spectral methods.

Inspired by [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl), it reimplements many of the same ideas with an emphasis on minimizing dependencies and providing a clean, composable interface for working with orthogonal polynomials, Fourier series, Hardy spaces, and other spectral bases.

## Features

- **Function approximation** using orthogonal polynomials (Jacobi, Ultraspherical, Hermite, Laguerre), Fourier series, Laurent series, and Hardy spaces
- **Operator algebra** — build differential, multiplication, and conversion operators symbolically, then materialize them as banded or dense matrices
- **Equation solving** — solve boundary value problems and operator equations using spectral collocation
- **Eigenvalue problems** — compute spectra of differential operators with high accuracy
- **Complex analysis** — Cauchy transforms, Hardy space projections, and Riemann–Hilbert problem solving
- **Adaptive methods** — automatic degree selection based on coefficient decay
- **Visualization** — built-in plotting of functions, coefficients, and domains

## Quick Example

Solve the Airy equation ``u'' - xu = 0`` on ``[-30, 30]`` with boundary conditions matching the Airy function:

```julia
using OperatorApproximation, Plots, SpecialFunctions

R = 30
sp  = Ultraspherical(0.0, UltraMappedInterval(-R, R, 2.0))
sp2 = Ultraspherical(2.0, UltraMappedInterval(-R, R, 2.0))

D  = Derivative()
M  = Multiplication(x -> x)
Op = D^2 - Conversion(sp2) * M

lbdry = FixedGridValues([-R], ChebyshevMappedInterval(-R, R)) |> Conversion
rbdry = FixedGridValues([ R], ChebyshevMappedInterval(-R, R)) |> Conversion

setbasis(sp)
u = (lbdry ⊘ rbdry ⊘ Op) \ [[airyai(-R)]; [airyai(R)]; x -> 0]
u = u[1]

plot(u; dx = 0.001)
```

## Installation

```julia
using Pkg
Pkg.add("OperatorApproximation")
```

Or install the development version directly:

```julia
using Pkg
Pkg.add(url = "https://github.com/tomtrogdon/OperatorApproximation.jl")
```

## Contents

```@contents
Pages = ["getting_started.md", "domains.md", "bases.md", "operators.md", "solvers.md", "examples.md", "api.md"]
Depth = 2
```

## Contributors

- Tom Trogdon (University of Washington)
- Kaitlynn Lilly
