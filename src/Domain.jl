abstract type Domain end
abstract type Interval <: Domain end
abstract type Circle <: Domain end
abstract type GridDomain end
abstract type GridRegion <: GridDomain end  # the grid is on the boundary

Tgrid = n -> cos.( (2*(1:n) .- 1)/(2*n) * pi ) |> reverse
Egrid = n -> cos.( (0:n-1)/(n-1) * pi ) |> reverse
LEgrid = n -> cos.( (1:n)/(n) * pi ) |> reverse
REgrid = n -> cos.( (0:n-1)/(n) * pi ) |> reverse

function DirectedEgrid(n)
    v = Egrid(n)
    v1 = ArgNum(v[1],0.0)
    vn = ArgNum(v[end],-pi)
    vcat([v1],v[2:end-1],[vn])
end

function DirectedLEgrid(n)
    v = LEgrid(n)
    v1 = ArgNum(v[1],0.0)
    vcat([v1],v[2:end])
end

function DirectedREgrid(n)
    v = REgrid(n)
    vn = ArgNum(v[end],-pi)
    vcat(v[1:end-1],[vn])
end

struct Exterior{T <: GridDomain} <: GridRegion
    D::Domain
    grid::Function
    GD::T
    function Exterior{T}(GD::T) where T <: GridDomain
        return new(GD.D,GD.grid,GD)
    end
end
Exterior(GD) = Exterior{typeof(GD)}(GD)

struct Interior{T <: GridDomain} <: GridRegion
    D::Domain
    grid::Function
    GD::T
    function Interior{T}(GD::T) where T <: GridDomain
        return new(GD.D,GD.grid,GD)
    end
end
Interior(GD) = Interior{typeof(GD)}(GD)


M = (A,B,x) -> (B - A)/2*x .+ (B + A)/2  # from I to [A,B]
iM = (A,B,x) -> 2/(B - A)*(x .- (B + A)/2) # From [A,B] to I

function isin(x::Number,I::Interval)
    X = I.imap(x)
    -1 <= real(X) <= 1 && imag(X) ≈ 0
end

function isin(x::Number,I::Circle)
    X = I.imap(x)
    abs2(X) ≈ 1
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

struct UnitCircle <: Circle  ## could be a map of UnitInterval()
    map::Function
    imap::Function
    cen::ComplexF64
    rad::Float64
    function UnitCircle()
        return new(x -> x, x -> x, 0.0im, 1.0)
    end
end

Base.show(io::IO, ::MIME"text/plain", z::UnitInterval)  =
           print(io, "UnitInterval(",z.a,",",z.b,")")

struct MappedInterval <: Interval
    map::Function # maps from I to interval
    imap::Function # the inverse map
    a
    b
    function MappedInterval(a,b)
        return new(x -> M(a,b,x), x -> iM(a,b,x), a, b)
    end
end

struct MappedCircle <: Circle
    map::Function # maps unit circle to mapped circle
    imap::Function
    cen::ComplexF64
    rad::Float64
    function MappedCircle(cen,rad)
        return new(x -> M(cen - rad, cen + rad,x), x -> iM(cen - rad, cen + rad,x), cen, rad)
    end
end

function ==(C1::Circle,C2::Circle)
    C1.cen ≈ C2.cen && C1.rad ≈ C2.rad
end

Base.show(io::IO, z::MappedInterval)  =
           print(io, "MappedInterval(",z.a,",",z.b,")")

function ==(I1::Interval,I2::Interval)
    I1.a ≈ I2.a && I1.b ≈ I2.b
end 

function iscompatible(GD1::GridDomain,GD2::GridDomain)
    GD1.D == GD2.D
end

abstract type GridInterval <: GridDomain end

struct ChebyshevInterval <: GridInterval
    D::Interval
    grid::Function
    function ChebyshevInterval()
        return new(UnitInterval(), Tgrid)
    end
end

struct LobattoInterval <: GridInterval
    D::Interval
    grid::Function
    function LobattoInterval()
        return new(UnitInterval(), Egrid)
    end
end

struct LLobattoInterval <: GridInterval # probably a bad name for this
    D::Interval
    grid::Function
    function LLobattoInterval()
        return new(UnitInterval(), Lgrid)
    end
end

struct RLobattoInterval <: GridInterval # probably a bad name for this
    D::Interval
    grid::Function
    function LLobattoInterval()
        return new(UnitInterval(), Rgrid)
    end
end

Base.show(io::IO, ::MIME"text/plain", z::ChebyshevInterval)  =
           print(io, "ChebyshevInterval(",sprint(print,z.D),")")

struct PeriodicInterval <: GridInterval
    D::Interval
    grid::Function
    function PeriodicInterval()
        return new(UnitInterval(), Pgrid)
    end
end

Base.show(io::IO, ::MIME"text/plain", z::PeriodicInterval)  =
           print(io, "PeriodicInterval(",sprint(print,z.D),")")

abstract type GridCircle <: GridDomain end 

struct PeriodicCircle <: GridCircle
    D::Circle
    grid::Function
    function PeriodicCircle()
        return new(UnitCircle(), Lgrid)
    end
end

struct JacobiInterval <: GridInterval
    D::Interval
    α::Number
    β::Number
    grid::Function
    function JacobiInterval(α,β)
        a, b = Jacobi_ab(α,β)
        gridfun = n -> Gauss_quad(a,b,n)[1]
        return new(UnitInterval(),α,β, gridfun)
    end
end

Base.show(io::IO, ::MIME"text/plain", z::JacobiInterval)  =
           print(io, "JacobiInterval(",sprint(print,z.D),",",z.α,",",z.β,")")

function ==(J1::JacobiInterval,J2::JacobiInterval)
    J1.D == J2.D && J1.α == J2.α && J1.β == J2.β
end

struct UltraInterval <: GridInterval
    D::Interval
    λ::Number
    grid::Function
    function UltraInterval(λ)
        a, b = Jacobi_ab(λ - 1/2,λ - 1/2)
        gridfun = n -> Gauss_quad(a,b,n-1)[1]
        return new(UnitInterval(), λ, gridfun)
    end
end

Base.show(io::IO, ::MIME"text/plain", z::UltraInterval)  =
           print(io, "UltraInterval(",sprint(print,z.D),",",z.λ,")")

function ==(J1::UltraInterval,J2::UltraInterval)
    J1.D == J2.D && J1.λ == J2.λ
end

struct ChebyshevMappedInterval <: GridInterval
    D::Interval
    grid::Function
    function ChebyshevMappedInterval(a,b)
        return new(MappedInterval(a,b), Tgrid)
    end
end

struct LobattoMappedInterval <: GridInterval
    D::Interval
    grid::Function
    function LobattoMappedInterval(a,b)
        return new(MappedInterval(a,b), Egrid)
    end
end

struct LLobattoMappedInterval <: GridInterval
    D::Interval
    grid::Function
    function LLobattoMappedInterval(a,b)
        return new(MappedInterval(a,b), Lgrid)
    end
end

struct RLobattoMappedInterval <: GridInterval
    D::Interval
    grid::Function
    function RLobattoMappedInterval(a,b)
        return new(MappedInterval(a,b), Rgrid)
    end
end

abstract type DirectedGridInterval <: GridInterval end

struct DirectedLobattoMappedInterval <: DirectedGridInterval
    D::Interval
    grid::Function
    dgrid::Function
    function DirectedLobattoMappedInterval(a,b)
        return new(MappedInterval(a,b), Egrid, DirectedEgrid)
    end
end

struct DirectedLLobattoMappedInterval <: DirectedGridInterval
    D::Interval
    grid::Function
    dgrid::Function
    function DirectedLLobattoMappedInterval(a,b)
        return new(MappedInterval(a,b), LEgrid, DirectedLEgrid)
    end
end

struct DirectedRLobattoMappedInterval <: DirectedGridInterval
    D::Interval
    grid::Function
    dgrid::Function
    function DirectedRLobattoMappedInterval(a,b)
        return new(MappedInterval(a,b), REgrid, DirectedREgrid)
    end
end


struct PeriodicMappedInterval <: GridInterval
    D::Interval
    grid::Function
    function PeriodicMappedInterval(a,b)
        return new(MappedInterval(a,b), Pgrid)
    end
end

struct PeriodicMappedCircle <: GridCircle
    D::Circle
    grid::Function
    function PeriodicMappedCircle(cen,rad)
        return new(MappedCircle(cen, rad), Lgrid)
    end
end

struct JacobiMappedInterval <: GridInterval
    D::Interval
    α::Float64
    β::Float64
    grid::Function
    function JacobiMappedInterval(a,b,α,β)
        A, B = Jacobi_ab(α,β)
        gridfun = n -> Gauss_quad(A,B,n-1)[1]
        return new(MappedInterval(a,b), α, β, gridfun)
    end
end

struct UltraMappedInterval <: GridInterval
    D::Interval
    λ::Float64
    grid::Function
    function UltraMappedInterval(a,b,λ)
        A, B = Jacobi_ab(λ - 1/2,λ - 1/2)
        gridfun = n -> Gauss_quad(A,B,n-1)[1]
        return new(MappedInterval(a,b), λ, gridfun)
    end
end