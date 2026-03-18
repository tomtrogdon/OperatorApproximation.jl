# Solving Equations

## Operator Equations

The `\` operator solves a system of the form `L * u = f` where `L` is a (block) `ConcreteOperator` and `f` is a right-hand side:

```julia
u = L \ f
```

The result is a vector of `BasisExpansion` objects (one per block row of the solution).

### Right-Hand Side Formats

The right-hand side can mix scalars, vectors, and anonymous functions:

```julia
u = (lbdry ⊘ rbdry ⊘ Op) \ [[val_left]; [val_right]; x -> rhs(x)]
```

Each entry corresponds to one block row:
- A scalar or length-1 vector satisfies the scalar boundary condition
- An anonymous function `x -> f(x)` is expanded in the range basis

## Boundary Value Problems

The standard pattern for a second-order BVP ``\mathcal{L}[u] = f`` on ``[a, b]`` with Dirichlet conditions ``u(a) = \alpha``, ``u(b) = \beta``:

```julia
sp  = Ultraspherical(0.0, ChebyshevMappedInterval(a, b))
sp2 = Ultraspherical(2.0, ChebyshevMappedInterval(a, b))

D   = Derivative()
Op  = D^2 + p * D + Conversion(sp2) * q   # example operator

lbdry = FixedGridValues([a], ChebyshevMappedInterval(a, b)) |> Conversion
rbdry = FixedGridValues([b], ChebyshevMappedInterval(a, b)) |> Conversion

setbasis(sp)
u = (lbdry ⊘ rbdry ⊘ Op) \ [[α]; [β]; x -> f(x)]
u = u[1]
```

The `⊘` operator stacks the rows: two boundary rows followed by the interior operator rows.

## Eigenvalue Problems

The `eigen` function computes the spectrum of a `ConcreteOperator` at a given truncation size `N`:

```julia
# Standard eigenvalue problem
vals, vecs = eigen(L, N)

# Generalized eigenvalue problem: L*u = λ M*u
vals, vecs = eigen(L, M, N)
```

The return type is `ContinuousEigen`, which packages eigenvalues with their corresponding `BasisExpansion` eigenvectors.

### Example: Harmonic Oscillator

```julia
sp  = HermiteFun(HermiteRealAxis())
sp2 = ...   # appropriate range space

D   = Derivative()
M   = Multiplication(x -> x^2)
H   = D^2 - Conversion(sp2) * M     # -∂² + x²

setbasis(sp)
ev = eigen(H, 100)
ev.values     # eigenvalues
ev.vectors    # eigenfunctions (BasisExpansion)
```

### Example: Hill's Equation

For the periodic eigenvalue problem (Floquet theory), use `FloquetDerivative`:

```julia
sp  = Fourier(PeriodicMappedInterval(0.0, 2π))
μ   = 0.3   # Floquet exponent

D   = FloquetDerivative(1, μ)
M   = Multiplication(x -> cos(2x))
Op  = D^2 + Conversion(sp2) * (λ - 2*cos(2x))

setbasis(sp)
ev = eigen(Op, 200)
```

## Riemann–Hilbert Problems

OperatorApproximation.jl includes a solver for Riemann–Hilbert problems (RHPs): find a matrix-valued function ``M(z)`` analytic off a contour ``\Gamma`` satisfying a jump condition ``M_+(z) = M_-(z) G(z)`` for a given jump matrix ``G``.

### Basic Interface

```julia
# Define the contour and jump matrix
D   = rhdomain(endpoints)
rng = rhrange(D)

# Build the RHP operator
L = rhmult(D, G) ⊘ rhrhs(D, G)

# Solve
sol = RHSolver(L, ...)
```

See the [Examples](@ref) page for complete worked examples.

## Tips for Convergence

- **Choose the right basis**: Ultraspherical/Chebyshev bases work well for smooth functions on compact intervals. Use Hermite/Laguerre for unbounded domains.
- **Match the range basis**: After applying `D^k`, the range is `Ultraspherical(λ + k, ...)`. Always use `Conversion(sp_range)` to bring operators back to a consistent basis.
- **Check `testconv`**: After solving, call `testconv(u)` to verify that the coefficient tail has decayed.
- **Increase `N`**: If convergence is poor, increase the truncation via `setN(N)` or by passing `N` explicitly to `eigen`.
