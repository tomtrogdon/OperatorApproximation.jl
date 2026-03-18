# Getting Started

## Installation

Install OperatorApproximation.jl from the Julia general registry:

```julia
using Pkg
Pkg.add("OperatorApproximation")
```

Then load it with:

```julia
using OperatorApproximation
```

## Core Concepts

OperatorApproximation.jl is built around three layered abstractions:

| Layer | Purpose |
|-------|---------|
| **Domain** | Describes where a function lives (an interval, circle, or the real axis) |
| **Basis** | Describes how functions are represented (Chebyshev, Ultraspherical, Fourier, …) |
| **Operator** | Describes an operation on functions (differentiation, multiplication, conversion, …) |

Functions are represented as [`BasisExpansion`](@ref) objects — a basis paired with a coefficient vector.

## Approximating a Function

Use the `BasisExpansion` constructor with a function and a basis. The package automatically determines the number of degrees of freedom needed for convergence:

```julia
using OperatorApproximation

# Chebyshev basis on [-1, 1]
sp = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))

# Approximate f(x) = sin(10x) adaptively
f = BasisExpansion(x -> sin(10x), sp)

# Evaluate at a point
f(0.5)

# Plot
plot(f)
```

You can also specify the number of terms explicitly:

```julia
f = BasisExpansion(x -> sin(10x), sp, 200)
```

## Setting a Default Basis

For convenience, set a global default basis so it need not be repeated:

```julia
setbasis(sp)

# Now BasisExpansion uses sp automatically
f = BasisExpansion(x -> exp(x))
```

## Differentiating a Function

Build a `Derivative` operator and apply it:

```julia
D  = Derivative()
sp = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))

setbasis(sp)
f  = BasisExpansion(x -> sin(x))
df = (D * f)[1]       # returns a BasisExpansion for cos(x)
```

## Multiplication by a Function

```julia
M  = Multiplication(x -> x)   # multiply by x
g  = (M * f)[1]               # x * sin(x)
```

## Solving a Boundary Value Problem

The operator `⊘` stacks operators vertically; `\` solves the resulting system:

```julia
sp  = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))
sp2 = Ultraspherical(2.0, ChebyshevMappedInterval(-1.0, 1.0))

D  = Derivative()
Op = D^2 + Conversion(sp2)      # u'' + u = f

# Boundary conditions: u(-1) = 0, u(1) = 0
lbdry = FixedGridValues([-1.0], ChebyshevMappedInterval(-1.0, 1.0)) |> Conversion
rbdry = FixedGridValues([ 1.0], ChebyshevMappedInterval(-1.0, 1.0)) |> Conversion

setbasis(sp)
u = (lbdry ⊘ rbdry ⊘ Op) \ [[0.0]; [0.0]; x -> cos(x)]
u = u[1]
```

## Convergence Testing

Check whether the coefficient vector of a `BasisExpansion` has converged:

```julia
testconv(f)    # prints convergence information
```

## Next Steps

- [Domains](@ref) — all supported domain types
- [Bases](@ref) — all supported basis types
- [Operators](@ref) — how to build and compose operators
- [Solving Equations](@ref) — boundary value problems, eigenvalue problems
- [Examples](@ref) — worked examples
