module OperatorApproximation

using SparseArrays, LinearAlgebra, Plots
import Plots: plot, +, -, *
import LinearAlgebra: I, Matrix

export Domain, GridDomain, Basis, Derivative, Evaluation, Ultraspherical, ChebyshevInterval, GridValues

include("Domain.jl")
include("Basis.jl")
include("Operators.jl")
include("Jacobi.jl")
include("Ultraspherical.jl")

end # module OperatorApproximation
