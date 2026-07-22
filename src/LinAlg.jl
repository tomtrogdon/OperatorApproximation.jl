# % LANCZOS Lanczos algorithm. --- https://www.cs.purdue.edu/archives/2002/wxg/codes/lanczos.m
# %
# %    Given the discrete inner product whose nodes are contained 
# %    in the first column, and whose weights are contained in the 
# %    second column, of the nx2 array xw, the call ab=LANCZOS(n,xw)
# %    generates the first n recurrence coefficients ab of the 
# %    corresponding discrete orthogonal polynomials. The n alpha-
# %    coefficients are stored in the first column, the n beta-
# %    coefficients in the second column, of the nx2 array ab.
# %
# %    The script is adapted from the routine RKPW in
# %    W.B. Gragg and W.J. Harrod, ``The numerically stable 
# %    reconstruction of Jacobi matrices from spectral data'', 
# %    Numer. Math. 44 (1984), 317-335.
# %
function lancz(N,x,w)
    Ncap = size(x,1);
    if N <= 0 || N > Ncap
        @error "lancz: N out of range"
    end
    p0 = copy(x)
    p1 = zeros(Ncap)
    p1[1] = w[1];
    for n = 1:Ncap-1
        pn = w[n+1]
        gam = 1
        sig=0
        t=0
        xlam = x[n+1];
        for k = 1:n+1
            rho = p1[k] + pn
            tmp = gam*rho
            tsig = sig
            if rho <= 0
                gam = 1
                sig = 0
            else
                gam = p1[k]/rho;
                sig = pn/rho;
            end
            tk = sig*(p0[k]-xlam) - gam*t
            p0[k] -= (tk-t)
            t = tk
            if sig <= 0
                pn = tsig*p1[k]
            else
            pn=(t^2)/sig;
            end
            tsig = sig
            p1[k] = tmp;
        end
    end
    SymTridiagonal(p0[1:N], sqrt.(p1[2:N]))
end

dot(f::BasisExpansion,g::BasisExpansion) = Base.sum(Multiplication(f)*Base.conj(g))

dot(f::BasisExpansion{T},g::BasisExpansion{T}) where T = LinearAlgebra.dot(f.c, g.c)

function dot(f::BasisExpansion{T},g::BasisExpansion{T}) where T <: DirectSum
    inner_prod = 0
    for i=1:length(f)
        inner_prod += dot(f[i],g[i])
    end
    return inner_prod
end

function sumdot(f::BasisExpansion,g::BasisExpansion)
    return dot(f,g)
end

#if f=[a,b], g=[c,d], this does <a,c> + <a,d> + <b,c> + <b,d>
function sumdot(f::BasisExpansion{T},g::BasisExpansion{T}) where T <: DirectSum
    inner_prod = 0
    for i=1:length(f)
        for j=1:length(g)
            inner_prod += dot(f[i],g[j])
        end
    end
    return inner_prod
end

function sumdot(f::BasisExpansion,g::BasisExpansion{T}) where T <: DirectSum
    inner_prod = 0
    for j=1:length(g)
        inner_prod += dot(f,g[j])
    end
    return inner_prod
end

function sumdot(f::BasisExpansion{T},g::BasisExpansion) where T <: DirectSum
    inner_prod = 0
    for i=1:length(f)
        inner_prod += dot(f[i],g)
    end
    return inner_prod
end

function sumdot(v1::Vector,v2::Vector)
    inner_prod = 0
    for j=1:length(v2)
        inner_prod += sumdot(v1[j],v2[j])
    end
    return inner_prod
end

# Block Lanczos + block Cholesky factorization.
#
# Mirrors the scalar SpikeDetect workflow
#   lanczos_with_extend → cholesky_jacobi → stieltjes_estimate / density
# but for p×p matrix orthogonal polynomials:
#   block_lanczos → block_cholesky_jacobi
#
# Block Lanczos on a Hermitian operator H with n×p starting block V₀
# produces a block Jacobi (block-tridiagonal) matrix
#
#   J_K = tridiag(α₁,…,αK ; β₁,…,βK₋₁)
#
# with p×p Hermitian diagonal blocks αk and p×p upper-triangular off-diagonal
# blocks βk (from QR).  The three-term recurrence is
#
#   H Vk = Vk₋₁ βk₋₁* + Vk αk + Vk₊₁ βk
#
# block_cholesky_jacobi then factors J_K = L L* where L is block lower-
# bidiagonal with diagonal blocks Γk (p×p lower triangular, from a p×p
# Cholesky) and subdiagonal blocks Δk (p×p general):
#
#   Γ₁ Γ₁* = α₁
#   Δk Γk* = βk          (solve for Δk given Γk and βk)
#   Γk Γk* = αk − Δk₋₁ Δk₋₁*   (k ≥ 2)
#
# The Cholesky factors (Γk, Δk) are the block analogue of the scalar (αk, βk)
# returned by cholesky_jacobi and are the natural input for the self-consistent
# density estimator (solve_mhat_analytic with CC = Γ∞ Γ∞*, CD = Δ∞ Γ∞*, DD = 0).

"""
    block_lanczos(apply_H, V0, K; reorth) -> (αs, βs)
    block_lanczos(H::AbstractMatrix, V0, K; reorth) -> (αs, βs)

Run K steps of block Lanczos on the Hermitian operator `H` (a matrix or a
function `v ↦ H*v`) with n×p starting block `V0`.

Returns
- `αs`: length-K vector of p×p Hermitian diagonal blocks
- `βs`: length-(K−1) vector of p×p upper-triangular off-diagonal blocks

The blocks satisfy the three-term recurrence

    H Vk = Vk₋₁ βk₋₁* + Vk αk + Vk₊₁ βk ,    βk upper triangular

and the block Jacobi matrix is J_K = tridiag(α₁,…,αK ; β₁,…,βK₋₁).

With `reorth=true` (default) double reorthogonalization against all previous
blocks is performed at each step for numerical stability.
"""
function block_lanczos(apply_H::Function, V0::AbstractMatrix, K::Integer;
                       reorth::Bool = true)
    n, p = size(V0)
    T = complex(float(eltype(V0)))

    Fqr0  = qr(Matrix{T}(V0))
    V     = Matrix(Fqr0.Q)[:, 1:p]          # orthonormalize starting block

    αs = Vector{Matrix{T}}(undef, K)
    βs = Vector{Matrix{T}}(undef, K - 1)
    Vs = Vector{Matrix{T}}(undef, K)
    Vs[1] = copy(V)

    Vprev = zeros(T, n, p)
    βprev = zeros(T, p, p)

    for k in 1:K
        HV = hcat(ntuple(j -> apply_H(@view V[:, j]), p)...)
        W  = HV - Vprev * βprev'

        α     = V' * W
        αs[k] = (α + α') / 2              # symmetrize numerically

        if k < K
            R = W - V * αs[k]

            if reorth
                for _ in 1:2, l in 1:k
                    R -= Vs[l] * (Vs[l]' * R)
                end
            end

            Fqr    = qr(R)
            Vnew   = Matrix(Fqr.Q)[:, 1:p]
            β      = Matrix(Fqr.R)[1:p, :]  # p×p upper triangular

            # enforce positive diagonal (canonical QR)
            for i in 1:p
                if real(β[i, i]) < 0
                    β[i, :]    .*= -1
                    Vnew[:, i] .*= -1
                end
            end

            βs[k]  = β
            Vprev  = V
            βprev  = β
            V      = Vnew
            if k < K - 1
                Vs[k + 1] = copy(V)
            end
        end
    end

    return αs, βs
end

function block_lanczos(H::AbstractMatrix, V0::AbstractMatrix, K::Integer; kwargs...)
    block_lanczos(v -> H * v, V0, K; kwargs...)
end

"""
    block_cholesky_jacobi(αs, βs) -> (Γs, Δs)

Block Cholesky factorization of the positive-definite block Jacobi matrix

    J_K = tridiag(α₁,…,αK ; β₁,…,βK₋₁)

produced by `block_lanczos`.  Returns (Γs, Δs) such that J_K = L L* where L
is block lower-bidiagonal with p×p lower-triangular diagonal blocks Γk and
p×p general subdiagonal blocks Δk.

Recurrence (analogous to the scalar `cholesky_jacobi`):

    Γ₁ Γ₁*  = α₁                       (p×p Cholesky)
    Δk  Γk* = βk                        (right-solve for Δk)
    Γk  Γk* = αk − Δk₋₁ Δk₋₁*          (k ≥ 2)

The factor L has the block structure

    L = [Γ₁          ]
        [Δ₁  Γ₂      ]
        [    Δ₂  Γ₃  ]
        [        ⋱  ⋱]

The Cholesky blocks (Γ∞, Δ∞) at convergence encode the bulk spectral
parameters: CC = Γ∞ Γ∞* and CD = Δ∞ Γ∞* match the coefficients of
`solve_mhat_analytic` with DD = 0.
"""
function block_cholesky_jacobi(αs::Vector, βs::Vector)
    K = length(αs)
    p = size(αs[1], 1)
    T = complex(float(eltype(αs[1])))

    Γs = Vector{Matrix{T}}(undef, K)
    Δs = Vector{Matrix{T}}(undef, K - 1)

    S = Matrix{T}(αs[1])                   # Schur complement  S = α₁ − Δ₀ Δ₀*  (Δ₀ = 0)
    Γs[1] = Matrix(cholesky(Hermitian(S)).L)

    for k in 1:K-1
        # Δk Γk* = βk  →  Δk = βk / Γk*  (right division by upper-triangular Γk*)
        Δs[k] = Matrix{T}(βs[k]) / Γs[k]'

        # Next Schur complement: αk₊₁ − Δk Δk*
        S      = Matrix{T}(αs[k+1]) - Δs[k] * Δs[k]'
        Γs[k+1] = Matrix(cholesky(Hermitian(S)).L)
    end

    return Γs, Δs
end