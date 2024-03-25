abstract type Basis end

abstract type CoefficientDomain end
struct ℤ <: CoefficientDomain end
struct ℕ₊ <: CoefficientDomain end
struct ℕ₋ <: CoefficientDomain end
struct 𝔼 <: CoefficientDomain end
struct 𝕏 <: CoefficientDomain end ## for when multiplication is not defined

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

getindex(b::DirectSum,i::Int64) = b.bases[i] |> DirectSum
getindex(b::DirectSum,i::UnitRange{Int64}) = b.bases[i] |> DirectSum
getindex(b::Basis,i) = getindex([b],i) |> DirectSum
axes(b::DirectSum) = axes(b.bases)
axes(b::DirectSum,i) = axes(b.bases,i)
axes(b::Basis) = axes([b])
axes(b::Basis,i) = axes([b],i)

function bases(b::DirectSum)
    b.bases
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
