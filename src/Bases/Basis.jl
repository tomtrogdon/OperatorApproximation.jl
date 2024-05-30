abstract type Basis end
abstract type FiniteBasis <: Basis end
abstract type CoefficientDomain end
struct ℤ <: CoefficientDomain end
struct ℕ₊ <: CoefficientDomain end
struct ℕ₋ <: CoefficientDomain end
struct 𝔼 <: CoefficientDomain end
struct 𝕏 <: CoefficientDomain end ## for when multiplication is not defined
struct AnyBasis <: Basis end

domainplot(b::Basis;kwargs...) = domainplot(b.GD.D;kwargs...)
domainplot!(b::Basis;kwargs...) = domainplot!(b.GD.D;kwargs...)
domainplot(V::Vector{T};kwargs...) where T <: Basis = domainplot([b.GD.D for b in V];kwargs...)

function ==(b1::AnyBasis,b2::Basis)
    true
end

function ==(b1::Basis,b2::AnyBasis)
    true
end

function ==(b1::AnyBasis,b2::AnyBasis)
    true
end

dim(b1::AnyBasis) = Inf

struct DirectSum <: Basis
    bases::Vector{T} where T <: Basis
    function DirectSum(b::Union{Vector{S},S}) where S <: Basis
        if !(typeof(b) <: Vector)
            return b
        elseif length(b) == 1
            return b[1]
        else
            new(b)
        end
    end
end
domainplot(b::DirectSum) = domainplot(b.bases)

getindex(b::DirectSum,i::Int64) = b.bases[i] |> DirectSum
getindex(b::DirectSum,i::UnitRange{Int64}) = b.bases[i] |> DirectSum
getindex(b::Basis,i) = getindex([b],i) |> DirectSum
axes(b::DirectSum) = axes(b.bases)
axes(b::DirectSum,i) = axes(b.bases,i)
axes(b::Basis) = axes([b])
axes(b::Basis,i) = axes([b],i)
length(b::DirectSum) = length(b.bases)


function bases(b::DirectSum)
    b.bases
end

function bases(b::Basis)
    [b]
end
function ⊕(b1::Basis,b2::Basis)
    DirectSum([b1,b2])
end

function ⊕(b1::Basis,b2::DirectSum)
    DirectSum(vcat([b1],b2.bases))
end

function ⊕(b1::DirectSum,b2::Basis)
    DirectSum(vcat(b1.bases,[b2]))
end

function ⊕(b1::DirectSum,b2::DirectSum)
    DirectSum(vcat(b1.bases,b2.bases))
end

function ⊕(bases...)
    if length(bases) == 1
        return bases[1]
    else
        out = bases[1] ⊕ bases[2]
        for i = 3:length(bases)
            out = out ⊕ bases[i]
        end
        return out
    end
end

function ==(b1::DirectSum,b2::DirectSum)
    prod(b1.bases .== b2.bases)
end

function isconvertible(b1::Basis,b2::Basis) # false by default
    false
end

#### SETTING UP A NEW BASIS ####
# (1) For each Basis, the isconvertible function should be 
#     overloaded to point to how it can be converted.
# (2) Each basis should have dim() implemented
# (3) If it makes sense, a routine to evaluate the
#     basis expansion should be implemented

include("BasisExpansion.jl")
include("GridValues.jl")
include("Jacobi.jl")
include("Ultraspherical.jl")
include("Fourier.jl")
include("Laurent.jl")
include("Hardy.jl")
include("Hermite.jl")
include("Erf.jl")
include("MarchenkoPastur.jl")

