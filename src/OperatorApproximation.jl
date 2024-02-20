module OperatorApproximation

using SparseArrays, LinearAlgebra, Plots
import Plots: plot
import Base: +, -, *, \, complex, /, length, iterate, log, sqrt, ==, ^
import LinearAlgebra: I, Matrix, norm

export Domain, GridDomain, Basis, Derivative, Evaluation, Ultraspherical, ChebyshevInterval,
     GridValues, FixedGridValues, FiniteGridValues, ConcreteOperator, Multiplication, ChebyshevMappedInterval, MappedInterval,
     LeftBoundaryFunctional, RightBoundaryFunctional, Projector, canonicalBC, BasisExpansion,
     BoundaryFunctional, Conversion, UnitInterval, MappedInterval, Transform, setbasis, setgrid, setN
     
struct StandardBasisVector
    j::Integer
end
    
function +(v::Vector,e::StandardBasisVector)
    w = copy(v)
    w[e.j] += 1
    return w
end

function +(e::StandardBasisVector,v::StandardBasisVector)
    w = copy(v)
    w[e.j] += 1
    return w
end
    
function -(v::Vector,e::StandardBasisVector)
     w = copy(v)
    w[e.j] -= 1
    return w
end

function -(e::StandardBasisVector,v::StandardBasisVector)
    w = copy(v)
    w[e.j] -= 1
    return w
end
    
function e(j)
    StandardBasisVector(j)
end
    
function e(j,n)
    fill(0.0,n) + e(j) 
end

## First handle domains, and approximation, transforms
include("Domain.jl")
include("Bases/Basis.jl")
include("Operators/AbstractOperators.jl")
include("Operators/ConcreteOperators.jl")
include("ArgNum.jl")
include("Cauchy.jl")
include("Solvers.jl")

global N = "adaptive"
function setN(n)
    global N = n
end

global basis = Ultraspherical(0.0,ChebyshevInterval())
function setbasis(b)
    global basis = b
end

global gridvals = GridValues(ChebyshevInterval())
function setgrid(b)
    global gridvals = b
end

Nmax = 10000
tol = 1e-14

end # module OperatorApproximation