module OperatorApproximation

using SparseArrays, LinearAlgebra, Plots
import Plots: plot, +, -, *
import LinearAlgebra: I, Matrix

export Domain, GridDomain, Basis, Derivative, Evaluation, Ultraspherical, ChebyshevInterval,
     GridValues, ConcreteOperator, Multiplication, ChebyshevMappedInterval, MappedInterval,
     LeftBoundaryFunctional, RightBoundaryFunctional

include("Domain.jl")
include("Basis.jl")
include("AbstractOperators.jl")
include("ConcreteOperators.jl")
include("Jacobi.jl")
include("Ultraspherical.jl")

end # module OperatorApproximation
