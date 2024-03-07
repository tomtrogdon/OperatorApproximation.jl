abstract type Basis end

struct DirectSum <: Basis
    bases::Vector{T} where T <: Basis
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
# (4) A transform should be implemented in the Operators 
#     directory corresponding to the basis.

include("BasisExpansion.jl")
include("GridValues.jl")
include("Jacobi.jl")
include("Ultraspherical.jl")
include("Fourier.jl")
