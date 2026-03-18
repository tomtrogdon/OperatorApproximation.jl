# Domains

A **domain** describes the underlying set on which functions are defined. Every basis in OperatorApproximation.jl is associated with a domain.

## Abstract Hierarchy

```
Domain
├── Interval          — a compact interval [a, b]
├── Circle            — a circle in the complex plane
├── Axis              — the full real line ℝ
├── SemiAxis          — a half-line [0, ∞)
└── DiscreteDomain    — a finite set of points (grid)
```

## Interval Domains

### Standard Intervals

| Constructor | Description |
|-------------|-------------|
| `ChebyshevMappedInterval(a, b)` | Interval `[a, b]` with Chebyshev (Gauss–Chebyshev) interior grid |
| `LobattoMappedInterval(a, b)` | Interval `[a, b]` with Chebyshev–Lobatto (boundary-inclusive) grid |
| `UltraMappedInterval(a, b, λ)` | Interval `[a, b]` with ultraspherical parameter `λ` |
| `JacobiMappedInterval(a, b, α, β)` | Interval `[a, b]` with Jacobi parameters `(α, β)` |

```julia
# Chebyshev interval on [-1, 1]
d = ChebyshevMappedInterval(-1.0, 1.0)

# Ultraspherical interval on [-30, 30], parameter λ = 2
d = UltraMappedInterval(-30.0, 30.0, 2.0)
```

### Periodic Intervals

| Constructor | Description |
|-------------|-------------|
| `PeriodicMappedInterval(a, b)` | Interval `[a, b]` with periodic (Fourier) grid |

## Circle Domains

| Constructor | Description |
|-------------|-------------|
| `UnitCircle()` | The unit circle `|z| = 1` |
| `MappedCircle(r)` | Circle of radius `r` |
| `PeriodicCircle()` | Unit circle with periodic parametrization |

Circle domains support interior (`Interior`) and exterior (`Exterior`) regions for complex analysis:

```julia
d = UnitCircle()
# Hardy{Interior, ...} — analytic inside the circle
# Hardy{Exterior, ...} — analytic outside the circle
```

## Axis and Semi-Axis Domains

| Constructor | Description |
|-------------|-------------|
| `HermiteRealAxis()` | Full real axis with Hermite weight `e^{-x²}` |
| `LaguerreSemiAxis()` | Half-axis `[0, ∞)` with Laguerre weight `e^{-x}` |
| `RationalRealAxis()` | Full real axis for rational function approximation |
| `RationalMappedAxis(a)` | Mapped real axis |

## Grid (Discrete) Domains

Discrete domains represent a finite set of interpolation or quadrature points. Most continuous domains have an associated grid domain obtained implicitly when constructing a `GridValues` basis.

| Constructor | Description |
|-------------|-------------|
| `FixedGridValues(pts, domain)` | Grid values at specified points `pts` |
| `GridValues` | Spectral grid values associated with a basis |

```julia
# Evaluate boundary condition at x = -1
lbdry = FixedGridValues([-1.0], ChebyshevMappedInterval(-1.0, 1.0)) |> Conversion
```

## Utility Functions

```julia
arclength(d)      # arc length of domain d
isin(x, d)        # true if point x lies in domain d
```
