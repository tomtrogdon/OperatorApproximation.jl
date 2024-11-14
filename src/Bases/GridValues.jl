abstract type DiscreteBasis <: Basis end # maybe change this name...

struct GridValues{T <: GridDomain} <: DiscreteBasis
    GD::T
end

struct FixedGridValues <: DiscreteBasis
    pts::Vector
    GD::GridDomain # Don't need the grid here, but for consistency...
    # All of the points should lie within GD
    # GD could be a continuous domain, or discrete
    function FixedGridValues(pts::Vector,GD::GridDomain)
        if  map(x -> isin(x,GD.D),pts) |> prod
            new(GD.D.imap.(pts),GD)
        else
            @error "Supplied points not in the Domain."
        end
    end
end
FixedGridValues(GD::Grid) = FixedGridValues(GD.grid,GD)

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
function dim(GV::FixedGridValues)
    GV.pts |> length
end

function dim(GV::GridValues)
    Inf
end

function pad(f::BasisExpansion{T},N) where T <: DiscreteBasis
    @error "Cannot pad a discrete basis"
    f
end

function testconv(f::BasisExpansion{T}) where T <: DiscreteBasis
    @warn "Cannot do a convergence test for a discrete basis"
    true
end

function chop(f::BasisExpansion{T}) where T <: DiscreteBasis
    @warn "Cannot chop discrete basis"
    f
end
####################################
#####  Important to implement  #####
####################################
function sum(f::BasisExpansion{T}) where T <: DiscreteBasis
    sum(f.c)
end

function moment(f::BasisExpansion{T},k::Int64) where T <: GridValues
    if k == 0
        return sum(f)
    end
    n = length(f.c)
    f.c.*(f.basis.GD.grid(n).^k) |> sum
end

function moment(f::BasisExpansion{T},k::Int64) where T <: FixedGridValues
    if k == 0
        return sum(f)
    end
    f.c.*(f.basis.pts.^k) |> sum
end


####################################
####################################
####################################

function (GV::GridValues)(N::Integer)
    GV.GD.D.map(GV.GD.grid(N))
    FixedGridValues(GV.GD.D.map(GV.GD.grid(N)),GV.GD)
end

function iscompatible(GV1::GridValues,GV2::GridValues)
    GV1 == GV2
end

function BasisExpansion(f::Function,basis::GridValues,N::Integer)
    grid = basis.GD.D.map(basis.GD.grid(N))
    BasisExpansion(basis,f.(grid))
end

#### FUNCTION OVERLOADING ####
function ^(f::BasisExpansion{T},a::Number) where T <: GridValues
    BasisExpansion(f.basis,(f.c).^a)
end

function ⊙(g::Function,f::BasisExpansion{T}) where T <: GridValues
    BasisExpansion(f.basis,g.(f.c))
end