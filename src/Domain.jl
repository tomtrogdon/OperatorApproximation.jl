abstract type Domain end
abstract type Interval <: Domain end
abstract type GridDomain end

M = (A,B,x) -> (B - A)/2*x .+ (B + A)/2  # from I to [A,B]
iM = (A,B,x) -> 2/(B - A)*(x .- (B + A)/2) # From [A,B] to I

function isin(x::Number,I::Interval)
    X = I.imap(x)
    -1 <= real(X) <= 1 && imag(X) ≈ 0
end

struct UnitInterval <: Interval
    map::Function
    imap::Function
    a
    b
    function UnitInterval()
        return new(x -> x, x -> x, -1.0, 1.0)
    end
end

struct MappedInterval <: Interval
    map::Function # maps from I to interval
    imap::Function # the inverse map
    a
    b
    function MappedInterval(a,b)
        return new(x -> M(a,b,x), x -> iM(a,b,x), a, b)
    end
end

function ==(I1::Interval,I2::Interval)
    I1.a ≈ I2.a && I1.b ≈ I2.b
end 

function iscompatible(GD1::GridDomain,GD2::GridDomain)
    GD1.D == GD2.D
end

struct ChebyshevInterval <: GridDomain
    D::Domain
    grid::Function
    function ChebyshevInterval()
        return new(UnitInterval(), Tgrid)
    end
end

struct JacobiInterval <: GridDomain
    D::Domain
    α::Number
    β::Number
    grid::Function
    function JacobiInterval(α,β)
        a, b = Jacobi_ab(α,β)
        gridfun = n -> Gauss_quad(a,b,n)
        return new(UnitInterval(),α,β, gridfun)
    end
end

function ==(J1::JacobiInterval,J2::JacobiInterval)
    J1.D == J2.D && J1.α == J2.α && J1.β == J2.β
end

struct UltraInterval <: GridDomain
    D::Domain
    λ::Number
    grid::Function
    function UltraInterval(λ)
        a, b = Jacobi_ab(λ - 1/2,λ - 1/2)
        gridfun = n -> Gauss_quad(a,b,n-1)[1]
        return new(UnitInterval(), λ, gridfun)
    end
end

function ==(J1::UltraInterval,J2::UltraInterval)
    J1.D == J2.D && J1.λ == J2.λ
end

struct ChebyshevMappedInterval <: GridDomain
    D::Domain
    grid::Function
    function ChebyshevMappedInterval(a,b)
        return new(MappedInterval(a,b), Tgrid)
    end
end

struct JacobiMappedInterval <: GridDomain
    D::Domain
    α::Float64
    β::Float64
    grid::Function
    function JacobiMappedInterval(a,b,α,β)
        A, B = Jacobi_ab(λ - 1/2,λ - 1/2)
        gridfun = n -> Gauss_quad(A,B,n-1)[1]
        return new(MappedInterval(a,b), α, β, gridfun)
    end
end

struct UltraMappedInterval <: GridDomain
    D::Domain
    λ::Float64
    grid::Function
    function UltraMappedInterval(a,b,λ)
        A, B = Jacobi_ab(λ - 1/2,λ - 1/2)
        gridfun = n -> Gauss_quad(A,B,n-1)[1]
        return new(MappedInterval(a,b), λ, gridfun)
    end
end