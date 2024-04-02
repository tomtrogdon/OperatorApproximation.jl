module OperatorApproximation

using SparseArrays, LinearAlgebra, Plots, FFTW, AbstractFFTs, HypergeometricFunctions
import Plots: plot, plot!
import Base: +, -, *, \, complex, /, length, iterate, log, sqrt, ==, ^,
    getindex, setindex!, firstindex, lastindex, show, getindex, size, axes,
    real, imag, abs
import LinearAlgebra: I, Matrix, norm, eigen, diagm

export Domain, GridDomain, Basis, Derivative, Evaluation, Ultraspherical, ChebyshevInterval,
     GridValues, FixedGridValues, FiniteGridValues, ConcreteOperator, Multiplication,
    ChebyshevMappedInterval, MappedInterval, BasisExpansion, Conversion, UnitInterval,
    MappedInterval, Transform, setbasis, setgrid, setN, UltraInterval, JacobiInterval,
    UltraMappedInterval, JacobiMappedInterval, PeriodicInterval, PeriodicMappedInterval,
    Fourier, Chop, eigen, plots, ⊞, DirectSum, ⊘, ⊕, FloquetDerivative, plot, plot!, eigen, \,
    ploteval, ploteval!, UnitCircle, MappedCircle, PeriodicCircle, PeriodicMappedCircle, Laurent, Hardy, Jacobi,
    CauchyTransform, Exterior, Interior, CauchyOperator, ArgNum, LobattoMappedInterval, LobattoInterval, BoundaryValue, BlockDiagonalAbstractOperator, AbstractZeroOperator, ZeroOperator,
    DirectedLobattoMappedInterval, DirectedLLobattoMappedInterval, DirectedRLobattoMappedInterval, Legendre, RHrange, RHdomain,
    BlockAbstractOperator, RHmult, RHrhs, matrix2BlockOperator

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
include("ArgNum.jl")
include("Operators/AbstractOperators.jl")
include("Operators/ConcreteOperators.jl")
include("Solvers.jl")
include("RHUtils.jl")

global N = "adaptive"
function setN(n)
    global N = n
end

global basis = Ultraspherical(0.0,ChebyshevInterval())
function setbasis(b)
    global basis = b
end

Nmax = 10000
tol = 1e-14

end # module OperatorApproximation