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

function getweight(sp::Ultraspherical)
    λ = sp.λ
    x -> JacobiW(λ - 0.5, λ - 0.5,x)
end
####################################
###### FUNCTION OVERLOADING ########
####################################
function sum(f::BasisExpansion{T}) where T <: Ultraspherical
    (f.basis.GD.D.b - f.basis.GD.D.a)*f.c[1]
end

function moment(f::BasisExpansion{T},k::Int64) where T <: Ultraspherical
    if k == 0
        return sum(f)
    end
    Multiplication(x -> x^k)*f |> sum
end

function (P::BasisExpansion{Ultraspherical})(X::Number) # Clenshaw's algorithm
    n = P.c |> length
    λ = P.basis.λ
    x = P.basis.GD.D.imap(X)
    if abs(imag(x)) > 1e-15 || abs(x) > 1
        return 0.0
    else
        a,b = Jacobi_ab(λ - 1/2, λ - 1/2)
        return (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
    end
end

function ^(f::BasisExpansion{T},x::Number) where T <: Ultraspherical
    # Should be done adaptively by padding if output coefs are not small
    a = f.basis.GD.D.a
    b = f.basis.GD.D.b
    gd = UltraMappedInterval(a,b,f.basis.λ)
    Conversion(f.basis)*((Conversion(GridValues(gd))*f)^x)
end

function ⊙(g::Function,f::BasisExpansion{T}) where T <: Ultraspherical
    # Should be done adaptively by padding if output coefs are not small
    a = f.basis.GD.D.a
    b = f.basis.GD.D.b
    gd = UltraMappedInterval(a,b,f.basis.λ)
    Conversion(f.basis)*(g ⊙ (Conversion(GridValues(gd))*f))
end

function norm(f::BasisExpansion{T}) where T <: Ultraspherical
    return norm(f.c)
end

# Comrade matrix
function roots(P::BasisExpansion{Ultraspherical}; tol = 1e-13)
    cs = chop(P.c)
    n = length(cs)
    n -= 1
    cs /= cs[end]
    λ = P.basis.λ
    a,b = Jacobi_ab(λ - 1/2, λ - 1/2)
    as = [a(i) for i in 0:n-1]
    bs = [b(i) for i in 0:n-1]  
    T = SymTridiagonal(as, bs[1:end-1]) |> Matrix
    T[end,1:end] -= cs[1:end-1]*bs[end]
    rtest = x -> abs(imag(x)) < tol && abs(real(x)) <= 1.0 + tol
    λ = eigvals(T)
    #λ[rtest.(λ)] |> real
    P.basis.GD.D.map.(λ[rtest.(λ)] |> real)
end