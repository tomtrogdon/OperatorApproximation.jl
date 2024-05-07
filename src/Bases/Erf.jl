#### Erf with decay at -inf ####
struct Erf <: FiniteBasis
    GD::GridAxis
end

cfd(sp::Erf) = ℕ₊

function dim(sp::Erf)
    1
end

function testconv(f::BasisExpansion{T}) where T <: Erf
    true
end

function chop(f::BasisExpansion{T}) where T <: Erf
    f
end

function (P::BasisExpansion{Erf})(X::Number) # Clenshaw's algorithm
    #x = P.basis.GD.D.imap(X)
    P.c[1]*0.5*(1 + erf(X/sqrt(2)))
end