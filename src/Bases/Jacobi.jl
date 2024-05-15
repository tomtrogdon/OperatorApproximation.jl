struct Jacobi <: Basis
    α::Number
    β::Number
    GD::GridInterval
end

function Legendre(a,b)
    GD = JacobiMappedInterval(a,b,0.0,0.0)
    Jacobi(0.0,0.0,GD)
end

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::Jacobi) = ℕ₊

function dim(sp::Jacobi)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Jacobi
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: Jacobi
    BasisExpansion(f.basis,chop(f.c))
end
####################################
#####  Important to implement  #####
####################################
function sum(f::BasisExpansion{T}) where T <: Jacobi
    (f.basis.GD.D.b - f.basis.GD.D.a)*f.c[1]
end

function moment(f::BasisExpansion{T},k::Int64) where T <: Jacobi
    if k == 0
        return sum(f)
    end
    Multiplication(x -> x^k)*f |> sum
end

function (P::BasisExpansion{Jacobi})(X::Number) # Clenshaw's algorithm
    n = P.c |> length
    α = P.basis.α
    β = P.basis.β
    x = P.basis.GD.D.imap(X)
    a,b = Jacobi_ab(α,β)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
end

