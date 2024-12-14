abstract type Domain end
abstract type Interval <: Domain end
abstract type Circle <: Domain end
abstract type Axis <: Domain end
abstract type SemiAxis <: Domain end
abstract type GridDomain end
abstract type GridRegion <: GridDomain end  # the grid is on the boundary

struct DiscreteDomain <: Domain
    map::Function
    imap::Function
    pts::Vector{ComplexF64}
    function DiscreteDomain(pts::Vector)
        new(x -> x, x -> x, pts)
    end
end
function isin(x::Number,D::DiscreteDomain)
    minimum(abs.( x .- D.pts)) < 1e-15
end

Tgrid = n -> cos.( (2*(1:n) .- 1)/(2*n) * pi ) |> reverse
Egrid = n -> cos.( (0:n-1)/(n-1) * pi ) |> reverse
LEgrid = n -> cos.( (1:n)/(n) * pi ) |> reverse
REgrid = n -> cos.( (0:n-1)/(n) * pi ) |> reverse

function DirectedEgrid(n)
    v = Egrid(n)
    v1 = ArgNum(v[1],1.0,0.0)
    vn = ArgNum(v[end],1.0,1.0*pi)
    vcat([v1],v[2:end-1],[vn])
end

function DirectedLEgrid(n)
    v = LEgrid(n)
    v1 = ArgNum(v[1],1.0,0.0)
    vcat([v1],v[2:end])
end

function DirectedREgrid(n)
    v = REgrid(n)
    vn = ArgNum(v[end],1.0,1.0*pi)
    vcat(v[1:end-1],[vn])
end

struct Exterior{T <: GridDomain} <: GridRegion
    D::Domain
    grid::Union{Function,Vector}
    GD::T
    function Exterior{T}(GD::T) where T <: GridDomain
        return new(GD.D,GD.grid,GD)
    end
end
Exterior(GD) = Exterior{typeof(GD)}(GD)

struct Interior{T <: GridDomain} <: GridRegion
    D::Domain
    grid::Union{Function,Vector}
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
    -1 <= real(X) <= 1 && imag(X) < 1e-14
end

function isin(x::Number,I::Circle)
    X = I.imap(x)
    abs(abs2(X) -1) < 1e-14
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

struct RealAxis <: Axis
    map::Function
    imap::Function
    cen::Union{ComplexF64,Float64}
    θ::Float64
    function RealAxis()
        return new(x -> x, x -> x, 0.0, 0.0)
    end
end

struct MappedAxis <: Axis
    map::Function
    imap::Function
    cen::Union{ComplexF64,Float64}
    θ::Float64
    function MappedAxis(σ,cen,θ)
        if θ ≈ 0.0
            return new(x -> (σ*x .+ cen), x -> (x .- cen)/σ, cen, θ)
        else
            return new(x -> (σ*x .+ cen)*exp(1im*θ), x -> (x*exp(-1im*θ).-cen)/σ, cen, θ)
        end
    end
end

struct PostiveRealAxis <: SemiAxis
    map::Function
    imap::Function
    cen::Union{ComplexF64,Float64}
    θ::Float64
    function PostiveRealAxis()
        return new(x -> x, x -> x, 0.0, 0.0)
    end
end

struct MappedSemiAxis <: SemiAxis
    map::Function
    imap::Function
    cen::Union{ComplexF64,Float64}
    θ::Float64
    function MappedSemiAxis(σ,cen,θ) # remove theta in favor of a complex number?
        # or stick with these special cases?
        if θ ≈ pi
            new(x -> -(σ*x .+ cen), x -> (-x.-cen)/σ, cen, θ)
        elseif θ ≈ 0.0
            new(x -> (σ*x .+ cen), x -> (x.-cen)/σ, cen, θ)
        else
            return new(x -> (σ*x .+ cen)*exp(1im*θ), x -> (x*exp(-1im*θ).-cen)/σ, cen, θ)
        end
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

arclength(I::Interval) = abs(I.a - I.b)
arclength(C::Circle) = 2*pi*C.rad

function ==(C1::Circle,C2::Circle)
    C1.cen ≈ C2.cen && C1.rad ≈ C2.rad
end

function ==(C1::Axis,C2::Axis)
    C1.cen ≈ C2.cen && C1.θ ≈ C2.θ
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
abstract type GridAxis <: GridDomain end
abstract type GridSemiAxis <: GridDomain end
abstract type GridCircle <: GridDomain end 

arclength(gd::GridDomain) = arclength(gd.D)

struct Grid <: GridDomain
    D::DiscreteDomain
    grid::Vector
end
Grid(D::DiscreteDomain) = Grid(D,D.pts)
Grid(V::Vector) = Grid(DiscreteDomain(V),V)

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

struct PeriodicCircle <: GridCircle
    D::Circle
    grid::Function
    function PeriodicCircle()
        return new(UnitCircle(), Lgrid)
    end
end

struct RationalRealAxis <: GridAxis
    D::Axis
    grid::Function
    function RationalRealAxis()
        Tm1 = z-> (1/1im)*((z.+1)./(z.-1)) #inverse mobius transform that maps unit circle onto real axis
        rat_mgrid = n-> (((0:n-1).+(1/2))./n).*(2*π) #[0,2π) shifted by 1/2 to avoid issues at 0 and infinity
        gridfun = n-> Tm1(exp.(1im.*rat_mgrid(n))) #x=T^{-1}(exp(iθ))
        return new(RealAxis(),gridfun)
    end
end

struct RationalMappedAxis <: GridAxis
    D::Axis
    grid::Function
    function RationalMappedAxis(σ,cen,θ)
        Tm1 = z-> (1/1im)*((z.+1)./(z.-1)) #inverse mobius transform that maps unit circle onto real axis
        rat_mgrid = n-> (((0:n-1).+(1/2))./n).*(2*π) #[0,2π) shifted by 1/2 to avoid issues at 0 and infinity
        gridfun = n-> Tm1(exp.(1im.*rat_mgrid(n))) #x=T^{-1}(exp(iθ))
        return new(MappedAxis(σ,cen,θ),gridfun)
    end
end

struct LaguerreSemiAxis <: GridSemiAxis
    D::SemiAxis
    grid::Function
    α::Float64
    function LaguerreSemiAxis(D,α)
        a, b = Laguerre_ab(α)
        gridfun = n -> Gauss_quad(a,b,n-1)[1]
        return new(D,gridfun,α)
    end
end

struct HermiteRealAxis <: GridAxis
    D::Axis
    grid::Function
    function HermiteRealAxis()
        a, b = Hermite_ab()
        gridfun = n -> Gauss_quad(a,b,n-1)[1]
        return new(RealAxis(),gridfun)
    end
end

struct HermiteAxis <: GridAxis
    ## TODO
end

struct JacobiInterval <: GridInterval
    D::Interval
    α::Number
    β::Number
    grid::Function
    function JacobiInterval(α,β)
        a, b = Jacobi_ab(α,β)
        gridfun = n -> Gauss_quad(a,b,n-1)[1]
        return new(UnitInterval(), α, β, gridfun)
    end
end

struct MarchenkoPasturInterval <: GridInterval
    D::Interval
    d::Number
    grid::Function
    function MarchenkoPasturInterval(d)
        a, b = MP_ab(d)
        gridfun = n -> Gauss_quad(a,b,n-1)[1]
        return new(MappedInterval((1 - sqrt(d))^2, (1 + sqrt(d))^2), d, gridfun)
    end
end

Base.show(io::IO, ::MIME"text/plain", z::JacobiInterval)  =
           print(io, "JacobiInterval(",sprint(print,z.D),",",z.α,",",z.β,")")

function ==(J1::JacobiInterval,J2::JacobiInterval)
    J1.D == J2.D && J1.α == J2.α && J1.β == J2.β
end

function ==(J1::MarchenkoPasturInterval,J2::MarchenkoPasturInterval)
    J1.D == J2.D && J1.d == J2.d
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

struct MarchenkoPasturMappedInterval <: GridInterval
    D::Interval
    d::Float64
    grid::Function
    function MarchenkoPasturMappedInterval(a,b,d)
        A, B = MP_ab(d)
        gridfun = n -> Gauss_quad(A,B,n-1)[1]
        return new(MappedInterval(a,b), d, gridfun)
    end
end
function ==(J1::MarchenkoPasturMappedInterval,J2::MarchenkoPasturMappedInterval)
    J1.D == J2.D && J1.d == J2.d
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