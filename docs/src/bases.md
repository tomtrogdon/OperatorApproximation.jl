# Bases

A **basis** associates a domain with a family of functions used to represent approximations. Each basis determines how coefficients relate to function values and which operators are available.

## BasisExpansion

All approximated functions are stored as `BasisExpansion` objects:

```julia
struct BasisExpansion{T<:Basis}
    basis :: T          # the basis
    c     :: Vector     # coefficient vector
end
```

### Constructors

```julia
# Adaptive: automatically choose N for convergence
f = BasisExpansion(x -> sin(x), sp)

# Fixed N
f = BasisExpansion(x -> sin(x), sp, 128)

# From coefficient vector
f = BasisExpansion(sp, c)
```

### Operations on BasisExpansion

| Operation | Description |
|-----------|-------------|
| `f(x)` | Evaluate at point `x` (or array of points) |
| `f + g`, `f - g` | Pointwise arithmetic (same basis) |
| `c * f` | Scalar multiplication |
| `sum(f)` | Definite integral (zeroth moment) |
| `moment(f, k)` | `k`-th moment |
| `norm(f)` | ℓ² norm of coefficients |
| `dot(f, g)` | Inner product |
| `pad(f, N)` | Truncate or zero-pad to `N` coefficients |
| `chop(f)` | Drop trailing small coefficients |
| `roots(f)` | Roots of the expansion |
| `testconv(f)` | Print convergence diagnostic |

## Polynomial Bases

### Ultraspherical

Ultraspherical (Gegenbauer) polynomials ``C_n^{(\lambda)}(x)`` with parameter ``\lambda``. The ``\lambda = 0`` case corresponds to Chebyshev polynomials of the first kind (up to normalization).

```julia
sp = Ultraspherical(λ, domain)

# Chebyshev-like basis on [-1, 1]
sp = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))

# Ultraspherical with λ=2 (common range basis after two derivatives)
sp2 = Ultraspherical(2.0, ChebyshevMappedInterval(-1.0, 1.0))
```

Ultraspherical bases are the workhorse for polynomial spectral methods. Taking one derivative maps ``C^{(\lambda)}`` to ``C^{(\lambda+1)}``, producing banded differentiation matrices.

### Jacobi

General Jacobi polynomials ``P_n^{(\alpha,\beta)}(x)`` with weight ``(1-x)^\alpha (1+x)^\beta``.

```julia
sp = Jacobi(α, β, domain)
sp = Jacobi(0.5, -0.5, ChebyshevMappedInterval(-1.0, 1.0))
```

### Legendre

Legendre polynomials (``\alpha = \beta = 0``):

```julia
sp = Jacobi(0.0, 0.0, ChebyshevMappedInterval(-1.0, 1.0))
```

## Fourier / Periodic Bases

### Fourier

Fourier series on a periodic interval. Coefficients are ordered as ``[a_0, a_1, b_1, a_2, b_2, \ldots]`` where ``a_k, b_k`` are cosine and sine amplitudes.

```julia
sp = Fourier(PeriodicMappedInterval(0.0, 2π))
f  = BasisExpansion(x -> exp(cos(x)), sp)
```

### Laurent

Laurent series on a circle in the complex plane. Natural for problems posed on the unit circle.

```julia
sp = Laurent(UnitCircle())
```

## Hardy Spaces

Hardy space bases for functions analytic inside or outside a disk. Used in complex analysis and Riemann–Hilbert problems.

```julia
# Functions analytic inside the unit circle
sp = Hardy{Interior}(UnitCircle())

# Functions analytic outside the unit circle
sp = Hardy{Exterior}(UnitCircle())
```

## Bases on Unbounded Domains

### Hermite

| Type | Description |
|------|-------------|
| `HermitePoly` | Hermite polynomials ``H_n(x)`` |
| `HermiteFun` | Hermite functions ``e^{-x^2/2} H_n(x)`` (square-integrable on ℝ) |

```julia
sp = HermiteFun(HermiteRealAxis())
f  = BasisExpansion(x -> exp(-x^2), sp)
```

### Laguerre

| Type | Description |
|------|-------------|
| `LaguerrePoly` | Laguerre polynomials ``L_n(x)`` on ``[0,\infty)`` |
| `LaguerreFun` | Laguerre functions ``e^{-x/2} L_n(x)`` |

```julia
sp = LaguerreFun(LaguerreSemiAxis())
```

### Oscillatory Rational (OscRational)

Oscillatory rational functions — useful for problems with oscillatory solutions on unbounded domains.

```julia
sp = OscRational(domain)
```

## Grid-Based Representations

`GridValues` stores a function by its values on a grid rather than by spectral coefficients. Useful as a building block for `Multiplication` operators and boundary conditions.

```julia
# Values at fixed points
lbdry = FixedGridValues([-1.0, 1.0], ChebyshevMappedInterval(-1.0, 1.0))
```

## Direct Sums

Bases can be combined into direct sums using `⊕` to handle piecewise or coupled systems:

```julia
sp1 = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 0.0))
sp2 = Ultraspherical(0.0, ChebyshevMappedInterval( 0.0, 1.0))
sp  = sp1 ⊕ sp2
```

## Setting a Default Basis

```julia
setbasis(sp)    # all subsequent BasisExpansion calls use sp by default
setN(N)         # set default number of degrees of freedom
```
