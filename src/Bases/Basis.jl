abstract type Basis end

function isconvertible(b1::Basis,b2::Basis) # false by default
    false
end

include("BasisExpansion.jl")
include("GridValues.jl")
include("Jacobi.jl")
include("Ultraspherical.jl")
