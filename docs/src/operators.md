# Operators

OperatorApproximation.jl uses a two-level operator design:

1. **`AbstractOperator`** — a symbolic description of an operation, independent of any particular basis.
2. **`ConcreteOperator`** — the matrix representation of an abstract operator acting between two specific bases.

This separation allows operators to be composed symbolically before any matrices are formed.

## Abstract Operators

### Derivative

```julia
D = Derivative()     # first derivative
D^2                  # second derivative (composition)
D^n                  # n-th derivative
```

For a `Ultraspherical(λ, ...)` basis, differentiation maps to `Ultraspherical(λ+1, ...)`:

```julia
sp  = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))
sp2 = Ultraspherical(1.0, ChebyshevMappedInterval(-1.0, 1.0))
# D maps sp → sp2 as a banded operator
```

### Multiplication

Multiply a function by another function:

```julia
M = Multiplication(x -> x)         # multiply by x
M = Multiplication(x -> sin(x))    # multiply by sin(x)
M = Multiplication(f)              # multiply by BasisExpansion f
```

For polynomial bases, multiplication by `x` is tridiagonal; multiplication by a general function is banded with bandwidth proportional to the expansion length of the multiplier.

### Conversion

Convert an expansion from one basis to another (when the bases are related by a change of basis):

```julia
Conversion(target_basis)
```

For Ultraspherical bases, `Conversion(Ultraspherical(λ+1, d))` is bidiagonal. Conversion operators are needed to reconcile the range of `Derivative` with a target basis:

```julia
sp  = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))
sp2 = Ultraspherical(2.0, ChebyshevMappedInterval(-1.0, 1.0))

D   = Derivative()
M   = Multiplication(x -> x)
Op  = D^2 - Conversion(sp2) * M    # u'' - x*u, maps sp → sp2
```

### Floquet Derivative

For problems with quasi-periodic solutions (Floquet theory), use:

```julia
FloquetDerivative(order, μ)    # ∂/∂x + iμ, order-th power
```

### Cauchy Transform

```julia
CauchyTransform()              # Cauchy singular integral operator
CauchyOperator(o)              # Cauchy operator with orientation o
```

### Boundary Values and Residues

```julia
BoundaryValue(o, range_basis)  # evaluation functional at boundary
Residue(range_basis)           # extract residues from poles
```

## Applying Operators

Applying an abstract operator to a basis produces a `ConcreteOperator`:

```julia
sp  = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))
setbasis(sp)

D   = Derivative()
cop = D * sp          # ConcreteOperator: sp → Ultraspherical(1.0, ...)
```

Applying a `ConcreteOperator` to a `BasisExpansion` returns a vector of `BasisExpansion` objects:

```julia
f  = BasisExpansion(x -> sin(x), sp)
df = (D * f)[1]       # first component: the derivative
```

## Composing Operators

### Arithmetic

```julia
Op1 + Op2      # sum of operators
Op1 - Op2      # difference
Op1 * Op2      # composition (Op2 applied first, then Op1)
c * Op         # scalar multiple
Op^n           # repeated composition
```

### Block Structure

Stack operators vertically with `⊘` and horizontally with `⊞`:

```julia
# Vertical stack (rows): boundary conditions above the differential operator
A = lbdry ⊘ rbdry ⊘ Op

# Horizontal concatenation (columns): side-by-side operators
B = Op1 ⊞ Op2

# Diagonal block operator
C = diagm([Op1, Op2, Op3])
```

### Direct Sum of Operators

For problems on direct-sum bases:

```julia
sp = sp1 ⊕ sp2
Op = matrix2BlockOperator([Op11 Op12; Op21 Op22])
```

## Truncation

Fix the size of an operator explicitly:

```julia
Truncation(Op, N)    # truncate operator to N×N
```

## Concrete Operators

A `ConcreteOperator` wraps a domain basis, range basis, and an underlying matrix (banded or dense). You can extract the matrix representation at a given size:

```julia
cop = D * sp
M   = Matrix(cop, N)    # N×N matrix
```

## Fast Operators

For large problems, fast variants avoid forming explicit matrices:

```julia
FastMultiplication(f)      # O(n log n) multiplication
FastConversion(target)     # fast basis conversion
```
