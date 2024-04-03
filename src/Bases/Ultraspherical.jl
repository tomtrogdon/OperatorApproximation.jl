struct Ultraspherical <: Basis
    λ::Number
    GD::GridInterval
end
####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::Ultraspherical) = ℕ₊

function dim(sp::Ultraspherical)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Ultraspherical
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: Ultraspherical  # add tolerance?
    BasisExpansion(f.basis,chop(f.c))
end
####################################
####################################
####################################

function (P::BasisExpansion{Ultraspherical})(X::Number) # Clenshaw's algorithm
    n = P.c |> length
    λ = P.basis.λ
    x = P.basis.GD.D.imap(X)
    a,b = Jacobi_ab(λ - 1/2, λ - 1/2)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
end