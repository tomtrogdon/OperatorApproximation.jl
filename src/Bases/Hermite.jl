abstract type Hermite <: Basis end

struct HermitePoly <: Hermite
    GD::GridAxis
end

struct HermiteFun <: Hermite
    GD::GridAxis
end

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::T) where T <: Hermite = ℕ₊

function dim(sp::T) where T <: Hermite
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Hermite
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: Hermite
    BasisExpansion(f.basis,chop(f.c))
end

function sum(f::BasisExpansion{T}) where T <: Hermite
    f.c[1]
end

function (P::BasisExpansion{HermitePoly})(X::Number) # Clenshaw's algorithm
    # probably should be weighted
    n = P.c |> length
    x = P.basis.GD.D.imap(X)
    a,b = Hermite_ab()
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
end

function (P::BasisExpansion{HermiteFun})(X::Number) # Clenshaw's algorithm
    # probably should be weighted
    n = P.c |> length
    x = P.basis.GD.D.imap(X)
    a,b = Hermite_ab()
    w = (2*pi)^(-.25)*exp(-x^2/4)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]*w
end