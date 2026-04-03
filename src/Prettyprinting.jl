########################################################################
#  PrettyPrinting.jl — human-readable show methods for
#  OperatorApproximation.jl
#
#  Convention:
#    show(io, x)                – compact, one-line form
#    show(io, MIME"text/plain", x) – multi-line / informational form
########################################################################

using Printf

# ─────────────────────────────────────────────────────────────────────
#  Helper utilities
# ─────────────────────────────────────────────────────────────────────

_fmt(x::Real)    = isinteger(x) ? string(Integer(x)) : string(x)
_fmt(x::Complex) = imag(x) == 0 ? _fmt(real(x)) :
                   real(x) == 0 ? string(imag(x))*"im" :
                   string(x)

"""Pretty-print an interval endpoint pair as [a, b]."""
_ivl(a, b) = "[" * _fmt(a) * ", " * _fmt(b) * "]"

"""Shorten a type name for display."""
_tname(x) = string(nameof(typeof(x)))

# ─────────────────────────────────────────────────────────────────────
#  1. Domains
# ─────────────────────────────────────────────────────────────────────

# --- Intervals ---

Base.show(io::IO, I::UnitInterval) =
    print(io, "UnitInterval ", _ivl(I.a, I.b))

Base.show(io::IO, ::MIME"text/plain", I::UnitInterval) =
    print(io, "UnitInterval ", _ivl(I.a, I.b))

Base.show(io::IO, I::MappedInterval) =
    print(io, "MappedInterval ", _ivl(I.a, I.b))

Base.show(io::IO, ::MIME"text/plain", I::MappedInterval) =
    print(io, "MappedInterval ", _ivl(I.a, I.b))

# --- Circles ---

Base.show(io::IO, C::UnitCircle) =
    print(io, "UnitCircle(center=", _fmt(C.cen), ", r=", _fmt(C.rad), ")")

Base.show(io::IO, ::MIME"text/plain", C::UnitCircle) =
    print(io, "UnitCircle  center = ", _fmt(C.cen), ", radius = ", _fmt(C.rad))

Base.show(io::IO, C::MappedCircle) =
    print(io, "MappedCircle(center=", _fmt(C.cen), ", r=", _fmt(C.rad), ")")

Base.show(io::IO, ::MIME"text/plain", C::MappedCircle) =
    print(io, "MappedCircle  center = ", _fmt(C.cen), ", radius = ", _fmt(C.rad))

# --- Axes ---

Base.show(io::IO, A::RealAxis) = print(io, "ℝ")

Base.show(io::IO, ::MIME"text/plain", A::RealAxis) =
    print(io, "RealAxis  ℝ  (center = ", _fmt(A.cen), ", θ = ", _fmt(A.θ), ")")

Base.show(io::IO, A::MappedAxis) =
    print(io, "MappedAxis(cen=", _fmt(A.cen), ", θ=", _fmt(A.θ), ")")

Base.show(io::IO, ::MIME"text/plain", A::MappedAxis) =
    print(io, "MappedAxis  center = ", _fmt(A.cen), ", θ = ", _fmt(A.θ))

# --- SemiAxes ---

Base.show(io::IO, A::PostiveRealAxis) = print(io, "ℝ₊")

Base.show(io::IO, ::MIME"text/plain", A::PostiveRealAxis) =
    print(io, "PostiveRealAxis  ℝ₊  (center = ", _fmt(A.cen), ", θ = ", _fmt(A.θ), ")")

Base.show(io::IO, A::MappedSemiAxis) =
    print(io, "MappedSemiAxis(cen=", _fmt(A.cen), ", θ=", _fmt(A.θ), ")")

Base.show(io::IO, ::MIME"text/plain", A::MappedSemiAxis) =
    print(io, "MappedSemiAxis  center = ", _fmt(A.cen), ", θ = ", _fmt(A.θ))

# --- DiscreteDomain ---

Base.show(io::IO, D::DiscreteDomain) =
    print(io, "DiscreteDomain(", length(D.pts), " pts)")

Base.show(io::IO, ::MIME"text/plain", D::DiscreteDomain) =
    print(io, "DiscreteDomain  ", length(D.pts), " points")

# ─────────────────────────────────────────────────────────────────────
#  2. Grid Domains
# ─────────────────────────────────────────────────────────────────────

# Helper to describe the underlying domain portion
_domain_str(gd::GridDomain) = sprint(show, gd.D)

# --- GridIntervals ---

Base.show(io::IO, gd::ChebyshevInterval) =
    print(io, "ChebyshevInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::ChebyshevInterval) =
    print(io, "ChebyshevInterval on ", gd.D)

Base.show(io::IO, gd::ChebyshevMappedInterval) =
    print(io, "ChebyshevMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::ChebyshevMappedInterval) =
    print(io, "ChebyshevMappedInterval on ", gd.D)

Base.show(io::IO, gd::LobattoInterval) =
    print(io, "LobattoInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::LobattoInterval) =
    print(io, "LobattoInterval on ", gd.D)

Base.show(io::IO, gd::LobattoMappedInterval) =
    print(io, "LobattoMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::LobattoMappedInterval) =
    print(io, "LobattoMappedInterval on ", gd.D)

Base.show(io::IO, gd::DirectedLobattoMappedInterval) =
    print(io, "DirectedLobattoMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::DirectedLobattoMappedInterval) =
    print(io, "DirectedLobattoMappedInterval on ", gd.D)

Base.show(io::IO, gd::DirectedLLobattoMappedInterval) =
    print(io, "DirectedLLobattoMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::DirectedLLobattoMappedInterval) =
    print(io, "DirectedLLobattoMappedInterval on ", gd.D)

Base.show(io::IO, gd::DirectedRLobattoMappedInterval) =
    print(io, "DirectedRLobattoMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::DirectedRLobattoMappedInterval) =
    print(io, "DirectedRLobattoMappedInterval on ", gd.D)

Base.show(io::IO, gd::PeriodicInterval) =
    print(io, "PeriodicInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::PeriodicInterval) =
    print(io, "PeriodicInterval on ", gd.D)

Base.show(io::IO, gd::PeriodicMappedInterval) =
    print(io, "PeriodicMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::PeriodicMappedInterval) =
    print(io, "PeriodicMappedInterval on ", gd.D)

Base.show(io::IO, gd::JacobiInterval) =
    print(io, "JacobiInterval(", gd.D, ", α=", _fmt(gd.α), ", β=", _fmt(gd.β), ")")

Base.show(io::IO, ::MIME"text/plain", gd::JacobiInterval) =
    print(io, "JacobiInterval on ", gd.D, "  α = ", _fmt(gd.α), ", β = ", _fmt(gd.β))

Base.show(io::IO, gd::JacobiMappedInterval) =
    print(io, "JacobiMappedInterval(", gd.D, ", α=", _fmt(gd.α), ", β=", _fmt(gd.β), ")")

Base.show(io::IO, ::MIME"text/plain", gd::JacobiMappedInterval) =
    print(io, "JacobiMappedInterval on ", gd.D, "  α = ", _fmt(gd.α), ", β = ", _fmt(gd.β))

Base.show(io::IO, gd::UltraInterval) =
    print(io, "UltraInterval(", gd.D, ", λ=", _fmt(gd.λ), ")")

Base.show(io::IO, ::MIME"text/plain", gd::UltraInterval) =
    print(io, "UltraInterval on ", gd.D, "  λ = ", _fmt(gd.λ))

Base.show(io::IO, gd::UltraMappedInterval) =
    print(io, "UltraMappedInterval(", gd.D, ", λ=", _fmt(gd.λ), ")")

Base.show(io::IO, ::MIME"text/plain", gd::UltraMappedInterval) =
    print(io, "UltraMappedInterval on ", gd.D, "  λ = ", _fmt(gd.λ))

Base.show(io::IO, gd::MarchenkoPasturInterval) =
    print(io, "MarchenkoPasturInterval(", gd.D, ", d=", _fmt(gd.d), ")")

Base.show(io::IO, ::MIME"text/plain", gd::MarchenkoPasturInterval) =
    print(io, "MarchenkoPasturInterval on ", gd.D, "  d = ", _fmt(gd.d))

Base.show(io::IO, gd::MarchenkoPasturMappedInterval) =
    print(io, "MarchenkoPasturMappedInterval(", gd.D, ", d=", _fmt(gd.d), ")")

Base.show(io::IO, ::MIME"text/plain", gd::MarchenkoPasturMappedInterval) =
    print(io, "MarchenkoPasturMappedInterval on ", gd.D, "  d = ", _fmt(gd.d))

# --- Grid on LLobattoMappedInterval / RLobattoMappedInterval ---

Base.show(io::IO, gd::LLobattoMappedInterval) =
    print(io, "LLobattoMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::LLobattoMappedInterval) =
    print(io, "LLobattoMappedInterval on ", gd.D)

Base.show(io::IO, gd::RLobattoMappedInterval) =
    print(io, "RLobattoMappedInterval(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::RLobattoMappedInterval) =
    print(io, "RLobattoMappedInterval on ", gd.D)

# --- GridCircles ---

Base.show(io::IO, gd::PeriodicCircle) =
    print(io, "PeriodicCircle(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::PeriodicCircle) =
    print(io, "PeriodicCircle on ", gd.D)

Base.show(io::IO, gd::PeriodicMappedCircle) =
    print(io, "PeriodicMappedCircle(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::PeriodicMappedCircle) =
    print(io, "PeriodicMappedCircle on ", gd.D)

# --- GridAxes ---

Base.show(io::IO, gd::RationalRealAxis) =
    print(io, "RationalRealAxis(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::RationalRealAxis) =
    print(io, "RationalRealAxis on ", gd.D)

Base.show(io::IO, gd::RationalMappedAxis) =
    print(io, "RationalMappedAxis(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::RationalMappedAxis) =
    print(io, "RationalMappedAxis on ", gd.D)

Base.show(io::IO, gd::HermiteRealAxis) =
    print(io, "HermiteRealAxis(", gd.D, ")")

Base.show(io::IO, ::MIME"text/plain", gd::HermiteRealAxis) =
    print(io, "HermiteRealAxis on ", gd.D)

# --- GridSemiAxes ---

Base.show(io::IO, gd::LaguerreSemiAxis) =
    print(io, "LaguerreSemiAxis(", gd.D, ", α=", _fmt(gd.α), ")")

Base.show(io::IO, ::MIME"text/plain", gd::LaguerreSemiAxis) =
    print(io, "LaguerreSemiAxis on ", gd.D, "  α = ", _fmt(gd.α))

# --- Grid (discrete) ---

Base.show(io::IO, gd::Grid) =
    print(io, "Grid(", length(gd.grid), " pts)")

Base.show(io::IO, ::MIME"text/plain", gd::Grid) =
    print(io, "Grid  ", length(gd.grid), " points on ", gd.D)

# --- Exterior / Interior ---

Base.show(io::IO, gd::Exterior) =
    print(io, "Exterior(", gd.GD, ")")

Base.show(io::IO, ::MIME"text/plain", gd::Exterior) =
    print(io, "Exterior of ", gd.GD)

Base.show(io::IO, gd::Interior) =
    print(io, "Interior(", gd.GD, ")")

Base.show(io::IO, ::MIME"text/plain", gd::Interior) =
    print(io, "Interior of ", gd.GD)

# ─────────────────────────────────────────────────────────────────────
#  3. Bases
# ─────────────────────────────────────────────────────────────────────

# --- Ultraspherical ---

Base.show(io::IO, sp::Ultraspherical) =
    print(io, "Ultraspherical(λ=", _fmt(sp.λ), ", ", sp.GD, ")")

function Base.show(io::IO, ::MIME"text/plain", sp::Ultraspherical)
    if sp.λ == 0.0
        name = "Chebyshev T"
    elseif sp.λ == 1.0
        name = "Chebyshev U"
    else
        name = "Ultraspherical(λ=" * _fmt(sp.λ) * ")"
    end
    print(io, name, " basis on ", sp.GD)
end

# --- Jacobi ---

Base.show(io::IO, sp::Jacobi) =
    print(io, "Jacobi(α=", _fmt(sp.α), ", β=", _fmt(sp.β), ", ", sp.GD, ")")

function Base.show(io::IO, ::MIME"text/plain", sp::Jacobi)
    if sp.α == 0.0 && sp.β == 0.0
        name = "Legendre"
    else
        name = "Jacobi(α=" * _fmt(sp.α) * ", β=" * _fmt(sp.β) * ")"
    end
    print(io, name, " basis on ", sp.GD)
end

# --- Fourier ---

Base.show(io::IO, sp::Fourier) =
    print(io, "Fourier(", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::Fourier) =
    print(io, "Fourier basis on ", sp.GD)

# --- Laurent ---

Base.show(io::IO, sp::Laurent) =
    print(io, "Laurent(", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::Laurent) =
    print(io, "Laurent basis on ", sp.GD)

# --- Hardy ---

function _hardy_label(::Type{T}) where T
    T <: Exterior ? "H⁻" :
    T <: Interior ? "H⁺" : "Hardy"
end

Base.show(io::IO, sp::Hardy{T,S}) where {T,S} =
    print(io, _hardy_label(T), "(", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::Hardy{T,S}) where {T,S} =
    print(io, _hardy_label(T), " Hardy space on ", sp.GD)

# --- Hermite ---

Base.show(io::IO, sp::HermitePoly) =
    print(io, "HermitePoly(", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::HermitePoly) =
    print(io, "Hermite polynomial basis on ", sp.GD)

Base.show(io::IO, sp::HermiteFun) =
    print(io, "HermiteFun(", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::HermiteFun) =
    print(io, "Hermite function basis on ", sp.GD)

# --- Laguerre ---

Base.show(io::IO, sp::LaguerrePoly) =
    print(io, "LaguerrePoly(α=", _fmt(sp.α), ", ", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::LaguerrePoly) =
    print(io, "Laguerre polynomial basis (α=", _fmt(sp.α), ") on ", sp.GD)

Base.show(io::IO, sp::LaguerreFun) =
    print(io, "LaguerreFun(α=", _fmt(sp.α), ", ", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::LaguerreFun) =
    print(io, "Laguerre function basis (α=", _fmt(sp.α), ") on ", sp.GD)

# --- Erf ---

Base.show(io::IO, sp::Erf) =
    print(io, "Erf(", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::Erf) =
    print(io, "Erf basis on ", sp.GD)

# --- OscRational ---

Base.show(io::IO, sp::OscRational) =
    print(io, "OscRational(α=", _fmt(sp.α), ", ", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::OscRational) =
    print(io, "OscRational basis (α=", _fmt(sp.α), ") on ", sp.GD)

# --- MarchenkoPastur ---

Base.show(io::IO, sp::MarchenkoPastur) =
    print(io, "MarchenkoPastur(d=", _fmt(sp.d), ", ", sp.GD, ")")

Base.show(io::IO, ::MIME"text/plain", sp::MarchenkoPastur) =
    print(io, "Marchenko–Pastur basis (d=", _fmt(sp.d), ") on ", sp.GD)

# --- Discrete Bases ---

Base.show(io::IO, gv::GridValues) =
    print(io, "GridValues(", gv.GD, ")")

Base.show(io::IO, ::MIME"text/plain", gv::GridValues) =
    print(io, "GridValues on ", gv.GD)

Base.show(io::IO, gv::FixedGridValues) =
    print(io, "FixedGridValues(", length(gv.pts), " pts, ", gv.GD, ")")

function Base.show(io::IO, ::MIME"text/plain", gv::FixedGridValues)
    print(io, "FixedGridValues  ", length(gv.pts), " point")
    length(gv.pts) != 1 && print(io, "s")
    print(io, " on ", gv.GD)
end

# --- AnyBasis ---

Base.show(io::IO, ::AnyBasis) = print(io, "AnyBasis")
Base.show(io::IO, ::MIME"text/plain", ::AnyBasis) = print(io, "AnyBasis (universal)")

# --- DirectSum of Bases ---

Base.show(io::IO, ds::DirectSum) =
    print(io, join([sprint(show, b) for b in ds.bases], " ⊕ "))

function Base.show(io::IO, ::MIME"text/plain", ds::DirectSum)
    println(io, "DirectSum of ", length(ds.bases), " bases:")
    for (i, b) in enumerate(ds.bases)
        prefix = i == length(ds.bases) ? "  └─ " : "  ├─ "
        i == length(ds.bases) ? print(io, prefix, sprint(show, MIME"text/plain"(), b)) :
                                println(io, prefix, sprint(show, MIME"text/plain"(), b))
    end
end

# ─────────────────────────────────────────────────────────────────────
#  4. BasisExpansion
# ─────────────────────────────────────────────────────────────────────

function Base.show(io::IO, f::BasisExpansion{T}) where T
    print(io, "BasisExpansion(", f.basis, ", ", length(f.c), " coeffs)")
end

function Base.show(io::IO, ::MIME"text/plain", f::BasisExpansion{T}) where T
    println(io, "BasisExpansion in ", sprint(show, MIME"text/plain"(), f.basis))
    if typeof(f.c[1]) <: Vector  # DirectSum case
        for (i, ci) in enumerate(f.c)
            i == length(f.c) ? print(io, "  block ", i, ": ", length(ci), " coefficients") :
                               println(io, "  block ", i, ": ", length(ci), " coefficients")
        end
    else
        n = length(f.c)
        println(io, "  ", n, " coefficient", n != 1 ? "s" : "")
        _show_coefficients(io, f.c)
    end
end

function _show_coefficients(io::IO, c::Vector; maxshow=6)
    n = length(c)
    if n <= maxshow
        for (i, v) in enumerate(c)
            i == n ? print(io, "    c[", i, "] = ", _fmt_coeff(v)) :
                     println(io, "    c[", i, "] = ", _fmt_coeff(v))
        end
    else
        half = maxshow ÷ 2
        for i in 1:half
            println(io, "    c[", i, "] = ", _fmt_coeff(c[i]))
        end
        println(io, "    ⋮")
        for i in n-half+1:n
            i == n ? print(io, "    c[", i, "] = ", _fmt_coeff(c[i])) :
                     println(io, "    c[", i, "] = ", _fmt_coeff(c[i]))
        end
    end
end

function _fmt_coeff(x::Number)
    if isa(x, Complex)
        r, im_ = reim(x)
        if abs(im_) < 1e-30
            return _fmt_real(r)
        elseif abs(r) < 1e-30
            return _fmt_real(im_) * "im"
        else
            return _fmt_real(r) * (im_ >= 0 ? " + " : " - ") * _fmt_real(abs(im_)) * "im"
        end
    else
        return _fmt_real(x)
    end
end

_fmt_real(x::Real) = abs(x) < 1e-30 ? "0" :
    isinteger(x) ? string(Integer(x)) :
    @sprintf("%.6g", x)

# ─────────────────────────────────────────────────────────────────────
#  5. Abstract Operators
# ─────────────────────────────────────────────────────────────────────

# --- Derivative ---

Base.show(io::IO, D::Derivative) =
    D.order == 1 ? print(io, "𝒟") : print(io, "𝒟", _superscript(D.order))

Base.show(io::IO, ::MIME"text/plain", D::Derivative) =
    D.order == 1 ? print(io, "Derivative (order 1)") :
    print(io, "Derivative (order ", D.order, ")")

# --- Shift ---

Base.show(io::IO, S::Shift) =
    S.order == 1 ? print(io, "𝒮") : print(io, "𝒮", _superscript(S.order))

Base.show(io::IO, ::MIME"text/plain", S::Shift) =
    S.order == 1 ? print(io, "Shift (order 1)") :
    print(io, "Shift (order ", S.order, ")")

# --- FloquetDerivative ---

Base.show(io::IO, D::FloquetDerivative) =
    print(io, "𝒟_F(μ=", _fmt(D.μ), ", order=", D.order, ")")

Base.show(io::IO, ::MIME"text/plain", D::FloquetDerivative) =
    print(io, "FloquetDerivative  order = ", D.order, ", μ = ", _fmt(D.μ))

# --- Multiplication ---

function _mul_label(f)
    if isa(f, Function)
        return "f(x)"
    elseif isa(f, BasisExpansion)
        return sprint(show, f)
    elseif isa(f, Vector)
        return "Vector(" * string(length(f)) * ")"
    else
        return string(f)
    end
end

Base.show(io::IO, M::Multiplication) =
    print(io, "M[", _mul_label(M.f), "]")

Base.show(io::IO, ::MIME"text/plain", M::Multiplication) =
    print(io, "Multiplication by ", _mul_label(M.f))

# --- Conversion / FastConversion / CoefConversion ---

Base.show(io::IO, C::Conversion) =
    print(io, "Conversion → ", C.range)

Base.show(io::IO, ::MIME"text/plain", C::Conversion) =
    print(io, "Conversion  → ", sprint(show, MIME"text/plain"(), C.range))

Base.show(io::IO, C::FastConversion) =
    print(io, "FastConversion → ", C.range)

Base.show(io::IO, ::MIME"text/plain", C::FastConversion) =
    print(io, "FastConversion  → ", sprint(show, MIME"text/plain"(), C.range))

Base.show(io::IO, C::CoefConversion) =
    print(io, "CoefConversion → ", C.range)

Base.show(io::IO, ::MIME"text/plain", C::CoefConversion) =
    print(io, "CoefConversion  → ", sprint(show, MIME"text/plain"(), C.range))

# --- BoundaryValue ---

Base.show(io::IO, B::BoundaryValue) =
    print(io, "BoundaryValue(order=", B.o, " → ", B.range, ")")

Base.show(io::IO, ::MIME"text/plain", B::BoundaryValue) =
    print(io, "BoundaryValue  order = ", B.o, " → ", sprint(show, MIME"text/plain"(), B.range))

# --- Residue ---

Base.show(io::IO, R::Residue) =
    print(io, "Residue → ", R.range)

Base.show(io::IO, ::MIME"text/plain", R::Residue) =
    print(io, "Residue  → ", sprint(show, MIME"text/plain"(), R.range))

# --- Truncation ---

Base.show(io::IO, T::Truncation) =
    print(io, "Truncation(", T.Op, ", k=", T.k, ")")

Base.show(io::IO, ::MIME"text/plain", T::Truncation) =
    print(io, "Truncation  k = ", T.k, " of ", sprint(show, MIME"text/plain"(), T.Op))

# --- CauchyTransform ---

Base.show(io::IO, ::CauchyTransform) = print(io, "CauchyTransform")
Base.show(io::IO, ::MIME"text/plain", ::CauchyTransform) = print(io, "Cauchy transform  𝒞")

# --- CauchyOperator ---

Base.show(io::IO, C::CauchyOperator) =
    print(io, "CauchyOperator(order=", C.o, ")")

Base.show(io::IO, ::MIME"text/plain", C::CauchyOperator) =
    C.o == 1 ? print(io, "Cauchy operator  𝒞⁺") :
    C.o == -1 ? print(io, "Cauchy operator  𝒞⁻") :
    print(io, "Cauchy operator  order = ", C.o)

# --- FourierTransform ---

Base.show(io::IO, F::FourierTransform) =
    print(io, "FourierTransform(order=", F.o, ")")

Base.show(io::IO, ::MIME"text/plain", F::FourierTransform) =
    print(io, "Fourier transform  order = ", F.o)

# --- AbstractZeroOperator ---

Base.show(io::IO, ::AbstractZeroOperator) = print(io, "𝟎")
Base.show(io::IO, ::MIME"text/plain", ::AbstractZeroOperator) = print(io, "Zero operator  𝟎")

# --- IdentityOperator ---

Base.show(io::IO, ::IdentityOperator) = print(io, "𝕀")
Base.show(io::IO, ::MIME"text/plain", ::IdentityOperator) = print(io, "Identity operator  𝕀")

# --- ProductOfAbstractOperators ---

function Base.show(io::IO, P::ProductOfAbstractOperators)
    join(io, [sprint(show, op) for op in P.Ops], " ∘ ")
end

function Base.show(io::IO, ::MIME"text/plain", P::ProductOfAbstractOperators)
    println(io, "Product of ", length(P.Ops), " operators:")
    for (i, op) in enumerate(P.Ops)
        prefix = i == length(P.Ops) ? "  └─ " : "  ├─ "
        i == length(P.Ops) ? print(io, prefix, sprint(show, op)) :
                             println(io, prefix, sprint(show, op))
    end
end

# --- SumOfAbstractOperators ---

function Base.show(io::IO, S::SumOfAbstractOperators)
    parts = String[]
    for (c, op) in zip(S.c, S.Ops)
        s = sprint(show, op)
        if c == 1
            push!(parts, s)
        elseif c == -1
            push!(parts, "-" * s)
        else
            push!(parts, string(c) * "⋅" * s)
        end
    end
    print(io, join(parts, " + "))
end

function Base.show(io::IO, ::MIME"text/plain", S::SumOfAbstractOperators)
    println(io, "Sum of ", length(S.Ops), " operators:")
    for (i, (c, op)) in enumerate(zip(S.c, S.Ops))
        prefix = i == length(S.Ops) ? "  └─ " : "  ├─ "
        coeff = c == 1 ? "" : c == -1 ? "-" : string(c) * " × "
        i == length(S.Ops) ? print(io, prefix, coeff, sprint(show, op)) :
                             println(io, prefix, coeff, sprint(show, op))
    end
end

# --- BlockAbstractOperator ---

function Base.show(io::IO, B::BlockAbstractOperator)
    n, m = size(B.Ops)
    print(io, n, "×", m, " BlockAbstractOperator")
end

function Base.show(io::IO, ::MIME"text/plain", B::BlockAbstractOperator)
    n, m = size(B.Ops)
    println(io, n, "×", m, " BlockAbstractOperator:")
    for i in 1:n
        row = ["  " * sprint(show, B.Ops[i, j]) for j in 1:m]
        i == n ? print(io, "  │ ", join(row, "  "), " │") :
                 println(io, "  │ ", join(row, "  "), " │")
    end
end

# --- BlockDiagonalAbstractOperator ---

function Base.show(io::IO, B::BlockDiagonalAbstractOperator)
    print(io, "diag(", join([sprint(show, op) for op in B.Ops], ", "), ")")
end

function Base.show(io::IO, ::MIME"text/plain", B::BlockDiagonalAbstractOperator)
    println(io, "BlockDiagonal of ", length(B.Ops), " operators:")
    for (i, op) in enumerate(B.Ops)
        prefix = i == length(B.Ops) ? "  └─ " : "  ├─ "
        i == length(B.Ops) ? print(io, prefix, sprint(show, op)) :
                             println(io, prefix, sprint(show, op))
    end
end

# ─────────────────────────────────────────────────────────────────────
#  6. Concrete Operators
# ─────────────────────────────────────────────────────────────────────

function Base.show(io::IO, Op::ConcreteOperator)
    print(io, "ConcreteOperator(", Op.domain, " → ", Op.range, ")")
end

function Base.show(io::IO, ::MIME"text/plain", Op::ConcreteOperator{D,R,T}) where {D,R,T}
    println(io, "ConcreteOperator")
    println(io, "  domain: ", sprint(show, MIME"text/plain"(), Op.domain))
    println(io, "  range:  ", sprint(show, MIME"text/plain"(), Op.range))
    print(io,   "  type:   ", _matrix_op_label(Op.L))
end

function _matrix_op_label(L::MatrixOperator)
    _tname(L)
end

function _matrix_op_label(L::SumOfMatrixOperators)
    "Sum of " * string(length(L.Ops)) * " matrix operators"
end

function _matrix_op_label(L::BlockMatrixOperator)
    n, m = size(L.Ops)
    string(n) * "×" * string(m) * " BlockMatrixOperator"
end

function _matrix_op_label(L::ZeroOperator)
    "ZeroOperator"
end

function _matrix_op_label(L::BasicBandedOperator)
    "BasicBandedOperator (bands " * string(L.nm) * ", " * string(L.np) * ")"
end

function _matrix_op_label(L::SemiLazyBandedOperator)
    "SemiLazyBandedOperator (bands " * string(L.nm) * ", " * string(L.np) * ")"
end

function _matrix_op_label(L::ProductOfBandedOperators)
    "Product of " * string(length(L.V)) * " banded operators"
end

function _matrix_op_label(L::TruncatedOperator)
    "TruncatedOperator (k=" * string(L.k) * ")"
end

function _matrix_op_label(L::FiniteRankOperator)
    "FiniteRankOperator (" * string(length(L.v)) * " rows)"
end

# ─────────────────────────────────────────────────────────────────────
#  7. ArgNum
# ─────────────────────────────────────────────────────────────────────

Base.show(io::IO, z::ArgNum) =
    print(io, "ArgNum(", _fmt(z.z), ", ρ=", _fmt(z.ρ), ", θ=", _fmt(z.θ), ")")

Base.show(io::IO, ::MIME"text/plain", z::ArgNum) =
    print(io, "ArgNum  z = ", _fmt(z.z), ", ρ = ", _fmt(z.ρ), ", θ = ", _fmt(z.θ))

# ─────────────────────────────────────────────────────────────────────
#  8. ContinuousEigen
# ─────────────────────────────────────────────────────────────────────

function Base.show(io::IO, E::ContinuousEigen)
    print(io, "ContinuousEigen(", length(E.values), " eigenvalues)")
end

function Base.show(io::IO, ::MIME"text/plain", E::ContinuousEigen)
    n = length(E.values)
    println(io, "ContinuousEigen  ", n, " eigenvalue", n != 1 ? "s" : "")
    maxshow = min(n, 6)
    for i in 1:maxshow
        λ = E.values[i]
        i == maxshow && n <= maxshow ? print(io, "  λ[", i, "] = ", _fmt_coeff(λ)) :
                                       println(io, "  λ[", i, "] = ", _fmt_coeff(λ))
    end
    n > maxshow && print(io, "  ⋮")
end

# ─────────────────────────────────────────────────────────────────────
#  9. RHP & RHSolver
# ─────────────────────────────────────────────────────────────────────

function Base.show(io::IO, rhp::RHP)
    n = size(rhp.Γ, 1)
    np = length(rhp.P)
    print(io, "RHP(", n, " contour", n != 1 ? "s" : "")
    np > 0 && print(io, ", ", np, " pole", np != 1 ? "s" : "")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", rhp::RHP)
    n = size(rhp.Γ, 1)
    np = length(rhp.P)
    nr = length(rhp.R)
    println(io, "Riemann–Hilbert Problem")
    println(io, "  ", n, " contour segment", n != 1 ? "s" : "")
    if np == 0 && nr == 0
        print(io, "  ", length(rhp.J), " jump matri", length(rhp.J) != 1 ? "ces" : "x")
    else
        println(io, "  ", length(rhp.J), " jump matri", length(rhp.J) != 1 ? "ces" : "x")
        if np > 0 && nr == 0
            print(io, "  ", np, " pole", np != 1 ? "s" : "")
        elseif np > 0
            println(io, "  ", np, " pole", np != 1 ? "s" : "")
        end
        nr > 0 && print(io, "  ", nr, " residue", nr != 1 ? "s" : "")
    end
end

function Base.show(io::IO, rhs::RHSolver)
    print(io, "RHSolver(", rhs.S.domain, " → ", rhs.S.range, ")")
end

function Base.show(io::IO, ::MIME"text/plain", rhs::RHSolver)
    nr = length(rhs.res)
    println(io, "RHSolver")
    println(io, "  operator: ", sprint(show, rhs.S))
    if nr == 0
        print(io, "  jumps:    ", length(rhs.jumps), " jump", length(rhs.jumps) != 1 ? "s" : "")
    else
        println(io, "  jumps:    ", length(rhs.jumps), " jump", length(rhs.jumps) != 1 ? "s" : "")
        print(io, "  residues: ", nr)
    end
end

# ─────────────────────────────────────────────────────────────────────
#  10. CoefficientDomain types
# ─────────────────────────────────────────────────────────────────────

Base.show(io::IO, ::ℤ) = print(io, "ℤ")
Base.show(io::IO, ::ℕ₊) = print(io, "ℕ₊")
Base.show(io::IO, ::ℕ₋) = print(io, "ℕ₋")
Base.show(io::IO, ::𝔼) = print(io, "𝔼")
Base.show(io::IO, ::𝕏) = print(io, "𝕏")

# ─────────────────────────────────────────────────────────────────────
#  Helper: integer superscripts for derivative / shift display
# ─────────────────────────────────────────────────────────────────────

const _SUPERSCRIPTS = Dict(
    '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
    '-' => '⁻')

function _superscript(n::Integer)
    s = string(n)
    String([get(_SUPERSCRIPTS, c, c) for c in s])
end