abstract type Laguerre <: Basis end

struct LaguerrePoly <: Laguerre
    GD::GridSemiAxis
    α::Float64
end

struct LaguerreFun <: Laguerre
    GD::GridSemiAxis
    α::Float64
end

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::T) where T <: Laguerre = ℕ₊

function dim(sp::T) where T <: Laguerre
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Laguerre
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: Laguerre
    BasisExpansion(f.basis,chop(f.c))
end

function sum(f::BasisExpansion{T}) where T <: Laguerre
    f.c[1]
end

function (P::BasisExpansion{LaguerrePoly})(X::Number) # Clenshaw's algorithm
    # probably should be weighted
    n = P.c |> length
    x = P.basis.GD.D.imap(X)
    α = P.basis.α
    a,b = Laguerre_ab(α)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
end

function (P::BasisExpansion{LaguerreFun})(X::Number) # Clenshaw's algorithm
    # probably should be weighted
    n = P.c |> length
    x = P.basis.GD.D.imap(X)
    α = P.basis.α
    a,b = Laguerre_ab(α)
    w = exp(-x/2)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]*w
end