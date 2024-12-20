module OperatorApproximation

using SparseArrays, LinearAlgebra, Plots, FFTW, AbstractFFTs, HypergeometricFunctions, Memoize, SpecialFunctions
import Plots: plot, plot!
import Base: +, -, *, \, complex, /, length, iterate, log, sqrt, ==, ^,
    getindex, setindex!, firstindex, lastindex, show, getindex, size, axes,
    real, imag, abs, sum, zero, intersect, conj, chop, copy
import LinearAlgebra: I, Matrix, norm, eigen, diagm, transpose, dot

export Domain, GridDomain, Basis, Derivative, Evaluation, Ultraspherical, ChebyshevInterval,
     GridValues, FixedGridValues, FiniteGridValues, ConcreteOperator, Multiplication, FastConversion,
    ChebyshevMappedInterval, MappedInterval, BasisExpansion, Conversion, UnitInterval,
    MappedInterval, Transform, setbasis, setgrid, setN, UltraInterval, JacobiInterval,
    UltraMappedInterval, JacobiMappedInterval, PeriodicInterval, PeriodicMappedInterval,
    Fourier, Chop, eigen, plots, ⊞, DirectSum, ⊘, ⊕, FloquetDerivative, plot, plot!, eigen, \,
    ploteval, ploteval!, UnitCircle, MappedCircle, PeriodicCircle, PeriodicMappedCircle, Laurent, Hardy, Jacobi,
    CauchyTransform, Exterior, Interior, CauchyOperator, ArgNum, LobattoMappedInterval, LobattoInterval, BoundaryValue, BlockDiagonalAbstractOperator, AbstractZeroOperator, ZeroOperator,
    DirectedLobattoMappedInterval, DirectedLLobattoMappedInterval, DirectedRLobattoMappedInterval, Legendre, rhrange, rhdomain,
    BlockAbstractOperator, rhmult, rhrhs, matrix2BlockOperator, RHSolver, domainplot, domainplot!, mvf2mof, coefplot, coefplot!, RHSolverVec,
    arclength, RHP, adapt, mofeval, mult2x2, dilog, rhwellposed, rhsplot, rhplot, clearCauchycache,
    HermitePoly, HermiteFun, Axis, GridAxis, RealAxis, HermiteRealAxis, CoefConversion, Erf, lancz, RecCoef,
    DiscreteDomain, Grid, Residue, moment, Truncation, MarchenkoPasturInterval, MarchenkoPastur, MarchenkoPasturMappedInterval,
    weightplot, weightplot!, RationalRealAxis, OscRational, dot, norm, ⊙, mfft, sumdot, combine, combinebasexp, simp,
    RationalMappedAxis, MappedSemiAxis, LaguerreSemiAxis, LaguerrePoly, LaguerreFun, roots,
    FourierTransform

function clearCauchycache()
    empty!(memoize_cache(cauchy))
    GC.gc()
end

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
include("OrthogonalPolynomials.jl")
include("Operators/AbstractOperators.jl")
include("Operators/ConcreteOperators.jl")
include("LinAlg.jl")
include("Solvers.jl")
include("RHUtils.jl")
include("Plotting.jl")


# probably not good to have these global variables
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