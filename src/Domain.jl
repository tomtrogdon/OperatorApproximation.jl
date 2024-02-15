abstract type Domain end
abstract type Interval <: Domain end
abstract type GridDomain end

M = (A,B,x) -> (B - A)/2*x .+ (B + A)/2  # from I to [A,B]
iM = (A,B,x) -> 2/(B - A)*(x .- (B + A)/2) # From [A,B] to I

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

struct ChebyshevInterval <: GridDomain
    D::Domain
    grid::Function
    function ChebyshevInterval()
        return new(UnitInterval(), Tgrid)
    end
end

struct ChebyshevMappedInterval <: GridDomain
    D::Domain
    grid::Function
    function ChebyshevMappedInterval(a,b)
        return new(MappedInterval(a,b), Tgrid)
    end
end