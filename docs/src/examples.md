# Examples

## 1. The Airy Equation

Solve the Airy equation ``u'' - xu = 0`` on ``[-30, 30]`` with boundary conditions from the exact solution:

```julia
using OperatorApproximation, Plots, SpecialFunctions

R = 30

# Domain and bases
sp  = Ultraspherical(0.0, UltraMappedInterval(-R, R, 2.0))
sp2 = Ultraspherical(2.0, UltraMappedInterval(-R, R, 2.0))

# Operators
D   = Derivative()
M   = Multiplication(x -> x)
Op  = D^2 - Conversion(sp2) * M

# Boundary conditions at x = ±R
lbdry = FixedGridValues([-R], ChebyshevMappedInterval(-R, R)) |> Conversion
rbdry = FixedGridValues([ R], ChebyshevMappedInterval(-R, R)) |> Conversion

# Solve
setbasis(sp)
u = (lbdry ⊘ rbdry ⊘ Op) \ [[airyai(-R)]; [airyai(R)]; x -> 0]
u = u[1]

plot(u; dx = 0.001)
```

The solution converges exponentially in the number of Ultraspherical coefficients. The mapped interval accounts for the oscillatory behavior near the turning point ``x = 0``.

---

## 2. Harmonic Oscillator Spectrum

Compute the eigenvalues of the quantum harmonic oscillator ``-u'' + x^2 u = \lambda u`` on the real line:

```julia
using OperatorApproximation

sp  = HermiteFun(HermiteRealAxis())
sp2 = ...   # range after two derivatives

D   = Derivative()
M   = Multiplication(x -> x^2)

setbasis(sp)
H  = D^2 - Conversion(sp2) * M
ev = eigen(H, 100)

# Exact eigenvalues: 2n+1 for n = 0, 1, 2, ...
ev.values[1:10]
```

---

## 3. Approximating an Analytic Function

Approximate ``f(z) = 1/(1 + 25z^2)`` (Runge's function) and plot its coefficient decay:

```julia
using OperatorApproximation, Plots

sp = Ultraspherical(0.0, ChebyshevMappedInterval(-1.0, 1.0))
f  = BasisExpansion(x -> 1 / (1 + 25x^2), sp)

# The coefficient plot shows geometric decay
coefplot(f)

# Evaluate at many points
xs = range(-1, 1; length = 1000)
ys = f.(xs)
plot(xs, ys)
```

---

## 4. Fourier Approximation of a Periodic Function

Approximate a periodic function and compute its derivative:

```julia
using OperatorApproximation, Plots

sp = Fourier(PeriodicMappedInterval(0.0, 2π))
f  = BasisExpansion(x -> exp(sin(x)), sp)

D  = Derivative()
setbasis(sp)
df = (D * f)[1]    # derivative: exp(sin(x)) * cos(x)

plot(f)
plot!(df)
```

---

## 5. Cauchy Transform on the Unit Circle

Compute the Cauchy transform of a function on the unit circle:

```julia
using OperatorApproximation

sp  = Laurent(UnitCircle())
f   = BasisExpansion(z -> 1 / (z - 2.0 + 0.0im), sp)

C   = CauchyTransform()
setbasis(sp)
Cf  = (C * f)[1]
```

---

## 6. Hill's Equation (Floquet Theory)

Compute the Floquet spectrum of the Mathieu equation ``u'' + (λ - 2q\cos(2t)) u = 0``:

```julia
using OperatorApproximation

q  = 0.5
sp = Fourier(PeriodicMappedInterval(0.0, 2π))

setbasis(sp)
# For a fixed Floquet exponent μ, use FloquetDerivative
μ   = 0.0
D   = FloquetDerivative(1, μ)
M   = Multiplication(t -> 2q * cos(2t))
# ... see the Hill's equation notebook for the complete setup
```

The notebooks directory contains a complete worked example for Hill's equation.

---

## 7. Nonlinear Schrödinger Equation via Riemann–Hilbert

The `notebooks/` directory contains a notebook demonstrating the solution of the defocusing NLS equation via the Riemann–Hilbert method implemented in OperatorApproximation.jl.

---

## Visualization

All `BasisExpansion` objects support direct plotting:

```julia
plot(f)                # line plot, auto dx
plot(f; dx = 0.001)    # specify point spacing

coefplot(f)            # log-scale coefficient magnitude

domainplot(d)          # visualize a domain

weightplot(sp)         # plot the basis weight function
```
