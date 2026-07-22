# Stochastic Lanczos Spectral Estimation
# Based on: "A Lanczos-Based Algorithmic Approach for Spike Detection
#            in Large Sample Covariance Matrices" (arXiv:2504.03066)
#
# Algorithms SR.1, SR.2, SR.3, P.2, P.3 from the paper.
# The (A ⊗ Iₖ) trick: run k independent Lanczos iterations on A with
# k random starting vectors (block structure).

"""
    lanczos_with_extend(apply_A, q1, n; tol, qwin) -> (a, b)

Run up to `n` steps of the Lanczos iteration with full reorthogonalization,
monitor convergence of the off-diagonal entries, and extend the output Jacobi
matrix with a constant tail once the entries have stabilised.

Convergence is declared when the standard deviation of the last `qwin`
off-diagonal entries falls below `tol * |mean| + tol`.  Early exit also occurs
if the off-diagonal norm drops below machine epsilon.

The returned arrays are extended by one entry: the final converged diagonal
value `a[end]` and off-diagonal value `b[end]` are each appended once more,
so the last two diagonal (and off-diagonal) entries are identical.  This
explicit repetition marks the semi-infinite constant tail used by
`stieltjes_estimate` and `density_components`.

Implements Algorithm B.1 of arXiv:2504.03066.
"""
function lanczos_with_extend(apply_A::Function, q1::AbstractVector, n::Integer;
                              tol::Real = 1e-6, qwin::Integer = 4)
    T = complex(eltype(q1))
    N = length(q1)
    a = zeros(real(T), n)
    b = zeros(real(T), n - 1)

    Q = Matrix{T}(undef, N, n)
    Q[:, 1] = convert(Vector{T}, q1) / norm(q1)

    j_conv = n   # step at which we declare convergence (default: full run)

    for j = 1:n
        Aq = apply_A(Q[:, j])
        a[j] = real(dot(Q[:, j], Aq))
        v = Aq - a[j] * Q[:, j]
        if j > 1
            v -= b[j-1] * Q[:, j-1]
        end

        # Full reorthogonalization against all previous vectors (twice for stability)
        for _ = 1:2
            for i = 1:j
                v -= dot(Q[:, i], v) * Q[:, i]
            end
        end

        β = norm(v)
        if j < n
            b[j] = β
        end

        # Early exit if Krylov space is exhausted
        if β < eps(real(T)) * n
            j_conv = j
            break
        end

        if j < n
            Q[:, j+1] = v / β
        end

        # Convergence check: std of last qwin off-diagonal entries is small
        if j >= qwin && j < n
            window = b[j-qwin+1:j]
            if std(window) < tol * abs(mean(window)) + tol
                j_conv = j
                break
            end
        end
    end

    # Extend: append one more copy of the converged tail values
    b_last = b[min(j_conv, n-1)]
    a_out  = vcat(a[1:j_conv], a[j_conv])
    b_out  = j_conv < n ? vcat(b[1:j_conv], b_last) : vcat(b, b_last, b_last)

    return a_out, b_out
end

"""
    cholesky_jacobi(a, b) -> (α, β)

Cholesky factorization J = L Lᵀ of a symmetric positive definite tridiagonal
(Jacobi) matrix with diagonal `a` and off-diagonal `b`.  Returns the
lower-bidiagonal Cholesky factor's diagonal `α` (length n) and sub-diagonal
`β` (length n-1).  Implements Algorithm B.2 of arXiv:2504.03066.
"""
function cholesky_jacobi(a::Vector, b::Vector)
    n = length(a)
    α = zeros(n)
    β = zeros(n - 1)
    α[1] = sqrt(max(a[1], 0.0))
    for i = 1:n-1
        α[i] == 0.0 && break
        β[i] = b[i] / α[i]
        α[i+1] = sqrt(max(a[i+1] - β[i]^2, 0.0))
    end
    return α, β
end

"""
    stieltjes_estimate(α, β, z) -> (m, γm, γp)

Compute the Stieltjes transform estimate m̂₀(z) and support endpoints (γ̂₋, γ̂₊)
from Cholesky factors `α` (length n) and `β` (length n-1) using the
continued-fraction backward recursion.  Implements Algorithm SR.2 of
arXiv:2504.03066.

The "tail" uses entries α[n-1] and β[n-1] (1-indexed), corresponding to the
constant part of the semi-infinite extension.
"""
function stieltjes_estimate(α::Vector{Float64}, β::Vector{Float64}, z::Number)
    n = length(α)   # β has length n-1
    n < 2 && error("Need at least 2 Cholesky entries")

    # Tail entries (0-indexed: n-2; 1-indexed: n-1)
    α_t = α[n-1]
    β_t = β[n-1]
    γm = (α_t - β_t)^2
    γp = (α_t + β_t)^2

    # Initialise Stieltjes transform at the tail (SR.2 step 2)
    m = (α_t^2 - z - β_t^2 +
         sqrt(complex(z - γp)) * sqrt(complex(z - γm))) / (2 * z * β_t^2)

    # Backward recursion from index n-2 to 1 (1-indexed); SR.2 step 3
    for i = n-2:-1:1
        αi = α[i]
        βi = β[i]
        m = 1 / (αi^2 - z - αi^2 * βi^2 * m / (1 + βi^2 * m))
    end

    return m, γm, γp
end

"""
    stieltjes_jacobi(a, b, z) -> (m, γm, γp)

Alternate Stieltjes transform estimator working directly from the Jacobi
(tridiagonal) matrix entries `a` (diagonal) and `b` (off-diagonal), as
returned by `lanczos_with_extend`.  Both vectors must have the same length n,
with the final two entries of each equal, marking the semi-infinite constant tail.

Support endpoints:  γ₋ = a∞ − 2b∞,  γ₊ = a∞ + 2b∞  (a∞ = a[end], b∞ = b[end]).

The Stieltjes transform m(z) = (J_ext − zI)⁻¹₁₁ is evaluated via the Schur
complement backward recursion (no Cholesky factorisation required):

    f_tail = ( (a∞ − z) + √((a∞ − z)² − 4b∞²) ) / (2b∞²)
    f_k    = 1 / (a[k] − z − b[k]² f_{k+1}),    k = n−1, …, 1
    m      = f₁
"""
function stieltjes_jacobi(a::Vector{Float64}, b::Vector{Float64}, z::Number)
    n = length(a)
    length(b) == n || throw(ArgumentError("a and b must have equal length (use lanczos_with_extend output)"))
    n ≥ 2 || error("need at least 2 entries")

    a∞ = a[n]
    b∞ = b[n]
    γm = a∞ - 2*b∞
    γp = a∞ + 2*b∞

    f = ((a∞ - z) - sqrt(complex((a∞ - z) - 2*b∞))*sqrt(complex((a∞ - z) + 2*b∞))) / (2*b∞^2)

    for k = n-1:-1:1
        f = 1 / (a[k] - z - b[k]^2 * f)
    end

    return f, γm, γp
end

"""
    cholesky_averaging!(all_α, all_β, navg)

In-place per-run tail stabilisation and cross-run tail synchronisation for a
collection of Cholesky factor pairs `(α, β)`.

For each pair, the last `navg` entries of the α and β tails are replaced by
their within-run means.  The per-run tail means are then averaged across all
runs, and every pair's tail region is overwritten with this common value so
all runs share the same asymptotic behaviour.

Implements Algorithm SR.3 of arXiv:2504.03066 as a standalone function.
"""
function cholesky_averaging!(all_α::Vector, all_β::Vector, navg::Integer)
    k = length(all_α)
    α_tails = zeros(k)
    β_tails = zeros(k)

    for j = 1:k
        αj = all_α[j];  nα = length(αj)
        βj = all_β[j];  nβ = length(βj)
        αt = mean(αj[max(1, nα - navg) : nα - 1])
        βt = mean(βj[max(1, nβ - navg + 1) : nβ])
        αj[max(1, nα - navg) : nα - 1] .= αt
        αj[nα] = αt
        βj[max(1, nβ - navg + 1) : nβ] .= βt
        α_tails[j] = αt
        β_tails[j] = βt
    end

    α_common = mean(α_tails)
    β_common = mean(β_tails)
    for j = 1:k
        αj = all_α[j];  nα = length(αj)
        βj = all_β[j];  nβ = length(βj)
        αj[max(1, nα - navg) : nα] .= α_common
        βj[max(1, nβ - navg + 1) : nβ] .= β_common
    end
    return nothing
end

function geomean(bvec)
    n = length(bvec)
    (prod(bvec))^(1/n)
end

"""
    jacobi_averaging!(all_a, all_b, navg)

In-place per-run tail stabilisation and cross-run tail synchronisation for a
collection of Jacobi matrix pairs `(a, b)` as returned by
`lanczos_with_extend` (both vectors have equal length, last two entries
identical).

Mirrors `cholesky_averaging!` but operates directly on the tridiagonal Jacobi
entries, so no Cholesky factorisation is required beforehand.
"""
function jacobi_averaging!(all_a::Vector, all_b::Vector, navg::Integer; sig = 0.01)
    k = length(all_a)
    a_tails = zeros(k)
    b_tails = zeros(k)

    for j = 1:k
        aj = all_a[j];  na = length(aj)
        bj = all_b[j];  nb = length(bj)
        at = mean(aj[max(1, na - navg) : na - 1])
        bt = geomean(bj[max(1, nb - navg) : nb - 1])
        aj[max(1, na - navg) : na - 1] .= at
        aj[na] = at
        bj[max(1, nb - navg) : nb - 1] .= bt
        bj[nb] = bt
        a_tails[j] = at
        b_tails[j] = bt
    end

    a_common = mean(a_tails)
    b_common = geomean(b_tails)*(1+sig) ## this is a hack.  Need to work it out.
    for j = 1:k
        aj = all_a[j];  na = length(aj)
        bj = all_b[j];  nb = length(bj)
        aj[max(1, na - navg) : na] .= a_common
        bj[max(1, nb - navg) : nb] .= b_common
    end
    return nothing
end

"""
    spike_detect(apply_A, N, k; nsteps, navg) -> (γm, γp, m0, cholesky_list)

Stochastic Lanczos spectral estimation with per-run tail stabilisation and
cross-run tail synchronisation.

For each of `k` independent Lanczos runs on the N-dimensional system:
1. Compute the Cholesky factorisation of the resulting Jacobi matrix.
2. Call `cholesky_averaging` to stabilise and synchronise the tails.

This avoids the near-degeneracy problem of the (A ⊗ Iₖ) trick, where
floating-point perturbations turn k exactly-repeated eigenvalues into k
nearly-equal but distinct ones.

# Arguments
- `apply_A`: function `ℝᴺ → ℝᴺ` representing the (symmetric) linear operator.
- `N`: ambient dimension.
- `k`: number of independent Lanczos runs.
- `nsteps`: Lanczos steps per run (default `⌈4 log N⌉`).
- `navg`: tail entries stabilised per run (default `⌈nsteps/4⌉`).

# Returns
- `γm`, `γp`: estimated left and right bulk-support edges (from run 1).
- `m0(z)`: Stieltjes transform estimate (from run 1).
- `cholesky_list`: length-`k` vector of `(α, β)` Cholesky factor pairs,
  all sharing a common stabilised tail. Use `cholesky_list[1]` for spike
  detection and density estimation.

Implements Algorithms SR.1–SR.3 and P.2 of arXiv:2504.03066.
"""
function spike_detect(apply_A::Function, N::Integer, k::Integer;
                      nsteps::Integer = max(5, ceil(Int, 4 * log(N))),
                      navg::Integer   = max(2, ceil(Int, nsteps / 4)))

    # ── SR.1: k independent Lanczos runs ─────────────────────────────────────
    all_α = Vector{Vector{Float64}}(undef, k)
    all_β = Vector{Vector{Float64}}(undef, k)
    for j = 1:k
        q = randn(N)
        q /= norm(q)
        a_jac, b_jac = lanczos_with_extend(apply_A, q, nsteps)
        all_α[j], all_β[j] = cholesky_jacobi(a_jac, b_jac)
    end

    # ── SR.3: tail stabilisation and cross-run synchronisation ───────────────
    cholesky_averaging!(all_α, all_β, navg)

    cholesky_list = [(all_α[j], all_β[j]) for j = 1:k]

    # ── SR.2: support endpoints and Stieltjes transform (from run 1) ─────────
    α1, β1 = cholesky_list[1]
    _, γm, γp = stieltjes_estimate(α1, β1, 1.0 + 1im)
    m0 = z -> stieltjes_estimate(α1, β1, z)[1]

    return γm, γp, m0, cholesky_list
end

"""
    count_spikes(γp, cholesky_list, N; C, δ) -> (count, outliers)

Count spike eigenvalues above the bulk edge using the first Cholesky pair in
`cholesky_list` (as returned by `spike_detect`).  Eigenvalues are obtained from
`eigvals(L'L)` where `L = Bidiagonal(α, β, :L)`, and those exceeding
`γp + C * N^(-δ)` are counted.

Implements Algorithm P.3 of arXiv:2504.03066.
"""
function count_spikes(γp::Float64,
                      cholesky_list::Vector{<:Tuple{Vector{Float64}, Vector{Float64}}},
                      N::Integer;
                      C::Float64 = 1.0, δ::Float64 = 0.25)
    threshold = γp + C * N^(-δ)
    α, β = cholesky_list[1]
    L = Bidiagonal(α, β, :L)
    λs = eigvals(L'*L)
    return count(λ -> λ > threshold, λs), λs[λs .> threshold]
end

support_endpoints(α_inf, β_inf) = ((α_inf - β_inf)^2, (α_inf + β_inf)^2)

"""
    density_components(α, β, x::Real) -> (u, v)

Real numbers (u, v) such that
    m₀(x + i0⁺) = u + i · v · √((b-x)(x-a))
for x in the bulk (a, b).

Each step of the continued-fraction recursion preserves this affine-in-S form:
if  m_{i+1} = u_{i+1} + i · v_{i+1} · S  with S² = (b-x)(x-a),  then setting

    P = 1 + β_i² · u_{i+1},     Q = β_i² · v_{i+1},     D = P² + Q²·S²,
    Ã = α_i² · P − x · D,        B̃ = α_i² · Q,
    denom = Ã² + B̃² · S²,

gives  u_i = Ã · D / denom,  v_i = B̃ · D / denom.

Both u and v are rational in x; this routine evaluates them with no complex
arithmetic. The result is also valid as the analytic continuation off the
bulk: for x ∉ [a,b], `denom` is the polynomial whose zeros locate the
discrete poles of m₀ (the spikes).
"""
function density_components(α::AbstractVector, β::AbstractVector, x::Real)
    N = length(α)
    @assert length(β) == N-1
    α∞, β∞ = α[N], β[N-1]
    a, b = support_endpoints(α∞, β∞)
    S2 = (b - x) * (x - a)

    # Base case: m_tail(x + i0⁺) = (α∞² − β∞² − x)/(2β∞²x) + i·(1/(2β∞²x))·S
    inv2βx = 1 / (2 * β∞^2 * x)
    u = (α∞^2 - β∞^2 - x) * inv2βx
    v = inv2βx

    @inbounds for i in (N-2):-1:1
        αi2 = α[i]^2
        βi2 = β[i]^2
        P = 1 + βi2 * u
        Q = βi2 * v
        D = P^2 + Q^2 * S2
        Ã = αi2 * P - x * D
        B̃ = αi2 * Q
        denom = Ã^2 + B̃^2 * S2
        u = Ã * D / denom
        v = B̃ * D / denom
    end
    return u, v
end

"""
    density(α, β, x::Real)

Density ρ(x) = v(x)·√((b−x)(x−a))/π on the bulk (a, b); zero outside.
"""
function density_est(α::AbstractVector, β::AbstractVector, x::Real)
    _, v = density_components(α, β, x)
    a, b = support_endpoints(α[end],β[end])
    S2 = (b - x) * (x - a)
    if S2 ≤ 0
        return zero(float(typeof(x)))
    end
    return v * sqrt(S2) / π
end

"""
    rational_factor(α, β, x::Real)

The rational factor r(x) such that ρ(x) = r(x)·√((b−x)(x−a)) on the bulk.
Defined for all real x; sign changes / blow-ups outside [a, b] mark the
discrete spike locations.
"""
function rational_factor(α::AbstractVector, β::AbstractVector, x::Real)
    _, v = density_components(α, β, x)
    return v / π
end

# Averaged overloads: accept the full cholesky_list from spike_detect and
# return the mean over all k Cholesky factors.
const CholeskyList = Vector{<:Tuple{Vector{Float64}, Vector{Float64}}}

function density_est(cholesky_list::CholeskyList, x::Real)
    mean(density_est(cj[1], cj[2], x) for cj in cholesky_list)
end

function rational_factor(cholesky_list::CholeskyList, x::Real)
    mean(rational_factor(cj[1], cj[2], x) for cj in cholesky_list)
end

"""
    solve_mhat_analytic(CC, CD, DD, z) -> (m, u, λ_inside)

Direct, non-iterative solution of the self-consistent equation

    m⁻¹ = CC - z·I - CD·(m⁻¹ + DD)⁻¹·CD*

on the branch with Im m ≻ 0. Coefficients: CC = ĈĈ*, CD = ĈD̂*, DD = D̂D̂*
(all Hermitian). Requires imag(z) > 0.

The substitution u = (m⁻¹ + DD)⁻¹ CD* yields the quadratic matrix equation

    CD·u² + (z·I - CC - DD)·u + CD* = 0,    m = (CC - z·I - CD·u)⁻¹.

Solved via the 2n×2n companion pencil. Because CD* = (CD)*, the 2n eigenvalues
come in reciprocal pairs (λ, 1/λ̄); for imag(z) > 0 exactly n lie inside the
unit disk. Branch selection takes the n eigenvalues of smallest modulus, which
is robust near the real axis where a geometric |λ|<1 threshold fails.
"""
function solve_mhat_analytic(CC::AbstractMatrix, CD::AbstractMatrix,
                             DD::AbstractMatrix, z::Number; disk_tol::Real=1e-9)
    imag(z) > 0 || throw(ArgumentError("need imag(z) > 0"))
    n  = LinearAlgebra.checksquare(CC)
    T  = complex(float(promote_type(eltype(CC), eltype(CD), eltype(DD), typeof(z))))
    In = Matrix{T}(I, n, n)
    Zn = zeros(T, n, n)

    CCt = Matrix{T}(CC); CDt = Matrix{T}(CD); DDt = Matrix{T}(DD)
    A  = CDt
    B  = z*In - CCt - DDt
    Cq = CDt'

    L0 = [Zn In; -Cq -B]
    L1 = [In Zn;  Zn  A]
    F  = eigen(L0, L1)
    λ  = F.values
    Vtop = F.vectors[1:n, :]

    p      = sortperm(abs.(λ))
    inside = p[1:n]
    gap    = abs(λ[p[n+1]]) - abs(λ[p[n]])
    gap > disk_tol || @warn("inside/outside moduli barely separated (gap=$gap): z is near the spectrum; use z + iη.")

    V = Vtop[:, inside]
    Λ = Diagonal(λ[inside])
    u = V * Λ * (V \ In)
    m = inv(CCt - z*In - CDt * u)
    return m, u, λ[inside]
end

"""
    bloch_matrix(Chat, Dhat, θ) -> Hermitian matrix

Symbol matrix 𝒜(θ) = (Ĉ - e^{-iθ} D̂)(Ĉ - e^{-iθ} D̂)*, Hermitian ≥ 0.
A real energy x lies in the a.c. spectrum iff x ∈ spec 𝒜(θ) for some θ ∈ [0, 2π).
"""
function bloch_matrix(Chat::AbstractMatrix, Dhat::AbstractMatrix, θ::Real)
    G = Chat .- cis(-θ) .* Dhat
    return Hermitian(G * G')
end

_bandvals(Chat, Dhat, θ) = eigvals(bloch_matrix(Chat, Dhat, θ))

function _vertex(_, x, k)
    (k == 1 || k == length(x)) && return x[k]
    d = x[k-1] - 2x[k] + x[k+1]
    abs(d) < 1e-14 && return x[k]
    return x[k] - (x[k-1] - x[k+1])^2 / (8d)
end

function _merge(iv; tol = 1e-7)
    isempty(iv) && return iv
    s = sort(iv; by = first)
    out = [collect(s[1])]
    for (a, b) in s[2:end]
        if a ≤ out[end][2] + tol
            out[end][2] = max(out[end][2], b)
        else
            push!(out, [a, b])
        end
    end
    return [(o[1], o[2]) for o in out]
end

"""
    spectral_support(Chat, Dhat; ngrid=4000)

Compute the support of the measure recovered from the self-consistent equation
via the band structure of the Bloch matrix 𝒜(θ), θ ∈ [0, 2π).

Returns a NamedTuple:
  intervals  :: Vector{Tuple}  — connected components of the support
  edges      :: Vector         — their endpoints (outer branch points)
  van_hove   :: Vector         — interior critical values (density singularities)
  θgrid, bands               — sampled grid and n sorted band curves (ngrid × n)
"""
function spectral_support(Chat::AbstractMatrix, Dhat::AbstractMatrix; ngrid::Integer = 4000)
    n = size(Chat, 1)
    θ = range(0, 2π; length = ngrid + 1)[1:end-1]
    W = Matrix{Float64}(undef, ngrid, n)
    for (i, t) in enumerate(θ)
        W[i, :] = _bandvals(Chat, Dhat, t)
    end

    intervals = Tuple{Float64,Float64}[]
    vanhove   = Float64[]
    for j in 1:n
        col = @view W[:, j]
        kmin = argmin(col); kmax = argmax(col)
        push!(intervals, (_vertex(θ, col, kmin), _vertex(θ, col, kmax)))
        d = diff(col)
        for k in 2:length(d)
            if sign(d[k-1]) * sign(d[k]) < 0
                push!(vanhove, _vertex(θ, col, k))
            end
        end
    end

    supp     = _merge(intervals)
    edges    = sort!(collect(Iterators.flatten(supp)))
    interior = sort!(filter(x -> any(a + 1e-6 < x < b - 1e-6 for (a, b) in supp), vanhove))

    return (intervals = supp, edges = edges, van_hove = interior,
            θgrid = collect(θ), bands = W)
end
