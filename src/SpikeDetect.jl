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
monitor convergence of the diagonal entries, and extend the output Jacobi
matrix with a constant tail once the entries have stabilised.

Convergence is declared when the standard deviation of the last `qwin`
diagonal entries falls below `tol * |mean| + tol`.  Early exit also occurs
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

        # Convergence check: std of last qwin diagonal entries is small
        if j >= qwin
            window = a[j-qwin+1:j]
            if std(window) < tol * abs(mean(window)) + tol
                j_conv = j
                break
            end
        end
    end

    # Extend: append one more copy of the converged tail values
    a_out = vcat(a[1:j_conv], a[j_conv])
    b_last = j_conv > 1 ? b[j_conv-1] : b[1]
    b_out  = vcat(b[1:j_conv-1], b_last, b_last)

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
    spike_detect(apply_A, N, k; nsteps, navg) -> (γm, γp, m0, jacobi)

Stochastic Lanczos spectral estimation using the (A ⊗ Iₖ) trick.

# Arguments
- `apply_A`: function `ℝᴺ → ℝᴺ` representing the (symmetric) linear operator.
- `N`: ambient dimension.
- `k`: number of blocks (columns of the starting matrix).
- `nsteps`: number of Lanczos iterations (default `⌈4 log N⌉`).
- `navg`: number of tail entries to average for tail stabilisation (SR.3,
  default `⌈nsteps/4⌉`).

# Starting-vector construction
Construct b₀ = [b₁; b₂; …; bₖ] ∈ ℝ^{Nk} where each bⱼ ∈ ℝᴺ is an
independent Gaussian vector normalised to unit 2-norm, then b₀ is normalised
to unit length.

# Matrix-vector product
`(A ⊗ Iₖ)` is applied implicitly: the current Nk vector is split into k
consecutive N-blocks and `apply_A` is called on each block independently.

# Returns
- `γm`, `γp`: estimated left and right edges of the bulk support.
- `m0(z)`: Stieltjes transform estimate (function of z ∈ ℂ).
- `jacobi`: 1-element vector containing the `(a, b)` Jacobi matrix pair from
  the single Lanczos run (for use with `count_spikes`).

Implements Algorithms SR.1–SR.3 and P.2 of arXiv:2504.03066.
"""
function spike_detect(apply_A::Function, N::Integer, k::Integer;
                      nsteps::Integer = max(5, ceil(Int, 4 * log(N))),
                      navg::Integer   = max(2, ceil(Int, nsteps / 4)))

    # ── Construct the Nk starting vector ────────────────────────────────────
    b0 = zeros(N * k)
    for j = 1:k
        bj = randn(N)
        bj /= norm(bj)                     # normalise each N-block to unit 2-norm
        b0[(j-1)*N+1 : j*N] = bj
    end
    b0 /= norm(b0)                         # normalise the full Nk vector

    # ── Implicit (A ⊗ Iₖ) matrix-vector product ─────────────────────────────
    function apply_kron(v::AbstractVector)
        w = similar(v)
        for j = 1:k
            idx = (j-1)*N+1 : j*N
            w[idx] = apply_A(v[idx])
        end
        return w
    end

    # ── SR.1: single Lanczos run on the Nk-dimensional system ───────────────
    a_jac, b_jac = lanczos_with_extend(apply_kron, b0, nsteps)
    jacobi = (a_jac, b_jac)

    # ── SR.3: stabilise the tail of the single Cholesky factorisation ───────
    α, β = cholesky_jacobi(a_jac, b_jac)
    nα = length(α)
    nβ = length(β)

    α_tail = mean(α[max(1, nα - navg) : nα - 1])
    β_tail = mean(β[max(1, nβ - navg + 1) : nβ])
    α[max(1, nα - navg) : nα - 1] .= α_tail
    β[max(1, nβ - navg + 1) : nβ]  .= β_tail
    α[nα] = α_tail   # enforce α[end] == α[end-1]
    cholesky = (α, β)

    # ── SR.2: Stieltjes transform and support endpoints ──────────────────────
    _, γm, γp = stieltjes_estimate(α, β, 1.0 + 1im)
    m0 = z -> stieltjes_estimate(α, β, z)[1]

    return γm, γp, m0, cholesky
end

"""
    count_spikes(γp, jacobi, N; C, δ) -> r̂

Count the number of spike eigenvalues above the bulk edge.  For each Lanczos
run, eigenvalues of the Jacobi matrix are computed via `eigvals(SymTridiagonal(a,b))`
and those exceeding the threshold `γp + C * N^(-δ)` are counted.  The returned
estimate is the rounded mean count across all k runs.

Implements Algorithm P.3 of arXiv:2504.03066.
"""
function count_spikes(γp::Float64,
                      cholesky::Tuple{Vector{Float64}, Vector{Float64}},
                      N::Integer;
                      C::Float64 = 1.0, δ::Float64 = 0.25)
    threshold = γp + C * N^(-δ)
    a, b = cholesky
    L = Bidiagonal(a,b,:L)
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
