# Domains

OperatorApproximation.jl separates the geometric description of a domain from its discretization into two distinct layers:

- A **`Domain`** describes the mathematical set (an interval, a circle, the real line, …). It carries only geometric information — endpoint values, maps, inverse maps — and knows nothing about how many points to use or where to put them.
- A **`GridDomain`** sits on top of a `Domain` by pairing it with a grid function `n -> Vector`. Given a size `n`, the grid function returns the `n` quadrature or interpolation points appropriate for a specific basis. All bases take a `GridDomain` as their argument, not a bare `Domain`.

```
Domain                          GridDomain
──────────────────────────      ──────────────────────────────────────────
Interval  ─────────────────>    GridInterval  (ChebyshevMappedInterval, …)
Circle    ─────────────────>    GridCircle    (PeriodicMappedCircle, …)
Axis      ─────────────────>    GridAxis      (HermiteRealAxis, …)
SemiAxis  ─────────────────>    GridSemiAxis  (LaguerreSemiAxis, …)
```

A `GridDomain` always stores:
- `.D` — the underlying `Domain`
- `.grid` — a function `n -> Vector` producing the grid of `n` points

---

## Domain Types

`Domain` types carry only geometric data. They are not passed directly to bases, but are embedded inside `GridDomain` instances.

### Interval Domains

| Type | Constructor | Description |
|------|-------------|-------------|
| `UnitInterval` | `UnitInterval()` | The reference interval ``[-1, 1]`` |
| `MappedInterval` | `MappedInterval(a, b)` | The interval ``[a, b]``, linearly mapped from ``[-1, 1]`` |
| `MappedCircle` | `MappedCircle(cen, rad)` | Circle of radius `rad` centered at `cen` |
| `UnitCircle` | `UnitCircle()` | Unit circle ``|z| = 1`` |
| `RealAxis` | `RealAxis()` | The full real line ``\mathbb{R}`` |
| `MappedAxis` | `MappedAxis(σ, cen, θ)` | Scaled/rotated real axis |
| `PostiveRealAxis` | `PostiveRealAxis()` | The half-line ``[0, \infty)`` |
| `MappedSemiAxis` | `MappedSemiAxis(σ, cen, θ)` | Scaled/rotated half-line |

Each domain stores a `map` (from the reference domain) and `imap` (its inverse), along with geometric data such as endpoints or center/radius.

---

## GridDomain Types

`GridDomain` types are what you actually pass to a basis constructor. Each one bundles a `Domain` with a specific grid rule.

### GridInterval

Discretizations of compact intervals:

| Constructor | Underlying grid | Typical use |
|-------------|-----------------|-------------|
| `ChebyshevMappedInterval(a, b)` | Chebyshev points of the first kind (interior) | Chebyshev / Ultraspherical(0) basis |
| `LobattoMappedInterval(a, b)` | Chebyshev–Lobatto points (includes endpoints) | Lobatto / Ultraspherical basis |
| `UltraMappedInterval(a, b, λ)` | Gauss–Ultraspherical quadrature for parameter `λ` | Ultraspherical(λ) basis |
| `JacobiMappedInterval(a, b, α, β)` | Gauss–Jacobi quadrature for parameters `(α, β)` | Jacobi(α, β) basis |
| `PeriodicMappedInterval(a, b)` | Uniform periodic grid | Fourier basis |
| `MarchenkoPasturMappedInterval(a, b, d)` | Gauss quadrature for Marchenko–Pastur weight | MarchenkoPastur basis |

```julia
# Chebyshev grid on [-1, 1] — suitable for Ultraspherical(0.0, ...) basis
gd = ChebyshevMappedInterval(-1.0, 1.0)
gd.D        # MappedInterval(-1.0, 1.0)  ← the geometry
gd.grid     # Tgrid  ← function: n -> n Chebyshev points

# Ultraspherical Gauss grid on [-30, 30] with λ = 2
gd = UltraMappedInterval(-30.0, 30.0, 2.0)
```

### Directed GridInterval

A subtype of `GridInterval` used in complex analysis contexts where orientation of the contour matters. Each directed grid tracks both a standard grid and a directed grid with tagged endpoints.

| Constructor | Description |
|-------------|-------------|
| `DirectedLobattoMappedInterval(a, b)` | Lobatto grid with both endpoints directed |
| `DirectedLLobattoMappedInterval(a, b)` | Lobatto grid with only left endpoint directed |
| `DirectedRLobattoMappedInterval(a, b)` | Lobatto grid with only right endpoint directed |

### GridCircle

Discretizations of circles in the complex plane:

| Constructor | Underlying grid | Typical use |
|-------------|-----------------|-------------|
| `PeriodicMappedCircle(cen, rad)` | Uniform periodic grid | Laurent series |

For complex analysis, `GridRegion` wraps a `GridDomain` with orientation information:

```julia
Exterior(gd)    # function is analytic outside the circle
Interior(gd)    # function is analytic inside the circle
```

These are used with `Hardy` bases.

### GridAxis

Discretizations of the real line:

| Constructor | Underlying grid | Typical use |
|-------------|-----------------|-------------|
| `HermiteRealAxis()` | Gauss–Hermite quadrature | Hermite polynomial/function basis |
| `RationalRealAxis()` | Möbius-transformed uniform periodic grid | Oscillatory rational basis |
| `RationalMappedAxis(σ, cen, θ)` | Scaled/rotated version of `RationalRealAxis` | Rational basis on a rotated axis |

### GridSemiAxis

Discretizations of semi-infinite domains:

| Constructor | Underlying grid | Typical use |
|-------------|-----------------|-------------|
| `LaguerreSemiAxis(D, α)` | Gauss–Laguerre quadrature with parameter `α` | Laguerre polynomial/function basis |

---

## Domain vs. GridDomain: Summary

| | `Domain` | `GridDomain` |
|---|---|---|
| Purpose | Geometric description | Geometric + discretization |
| Fields | `map`, `imap`, endpoints | `.D` (Domain), `.grid` (Function) |
| Passed to bases? | No | Yes |
| Examples | `MappedInterval(a,b)` | `ChebyshevMappedInterval(a,b)` |

---

## Utility Functions

```julia
arclength(d)      # arc length of a Domain or GridDomain
isin(x, d)        # true if point x lies in domain d
iscompatible(gd1, gd2)   # true if gd1 and gd2 share the same underlying Domain
```
