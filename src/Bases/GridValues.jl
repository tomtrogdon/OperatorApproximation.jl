abstract type DiscreteBasis <: Basis end # maybe change this name...

struct GridValues <: DiscreteBasis
    GD::GridDomain
end

function dim(GV::GridValues)
    Inf
end

# struct FiniteGridValues <: DiscreteBasis  # is this needed?
#     N::Integer
#     GD::GridDomain
# end

# function dim(GV::FiniteGridValues)
#     GV.N
# end

struct FixedGridValues <: DiscreteBasis
    pts::Vector
    GD::GridDomain # Don't need the grid here, but for consistency...
    function FixedGridValues(pts::Vector,GD::GridDomain)
        if  map(x -> isin(x,GD.D),pts) |> prod
            new(GD.D.imap.(pts),GD)
        else
            @error "Supplied points not in the Domain."
        end
    end
end

function (GV::GridValues)(N::Integer)
    GV.GD.D.map(GV.GD.grid(N))
    FixedGridValues(GV.GD.D.map(GV.GD.grid(N)),GV.GD)
end

function dim(GV::FixedGridValues)
    GV.pts |> length
end

function iscompatible(GV1::GridValues,GV2::GridValues)
    GV1 == GV2
end

function BasisExpansion(f::Function,basis::GridValues,N::Integer)
    grid = basis.GD.D.map(basis.GD.grid(N))
    BasisExpansion(basis,f.(grid))
end