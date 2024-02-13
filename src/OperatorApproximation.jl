module OperatorApproximation

using SparseArrays, LinearAlgebra, Plots
import Plots: plot, plot!
import Base: +, -, *, \, length, sqrt, /, iterate, log
import LinearAlgebra: I, Matrix, norm

export Domain, GridDomain, Basis, Derivative, Evaluation, Ultraspherical, ChebyshevInterval,
     GridValues, ConcreteOperator, Multiplication, ChebyshevMappedInterval, MappedInterval,
     LeftBoundaryFunctional, RightBoundaryFunctional, Projector, canonicalBC, BasisExpansion,
     BoundaryFunctional, setN, setbasis

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

include("Domain.jl")
include("Basis.jl")
include("ArgNum.jl")
include("AbstractOperators.jl")
include("ConcreteOperators.jl")
include("BasisExpansion.jl")
include("Jacobi.jl")
include("Cauchy.jl")
include("Ultraspherical.jl")

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

function \(L::Vector{ConcreteOperator},b::Vector,N::Integer)
    Ops = []; rhss = []
    for i = 1:length(L)
        if typeof(L[i]) <: ConcreteFunctional
            push!(Ops,Matrix(L[i],N))
            push!(rhss,b[i])
        end
    end
    Op = vcat(Ops...)
    rhs = vcat(rhss...)
    k = size(Op)[1]
    ## Assume only one non Functional Operator
    for i = 1:length(L)
        if !(typeof(L[i]) <: ConcreteFunctional)
            Op = vcat(Op,Matrix(L[i],N-k,N))
            P = Projector(N-k)
            rhs = vcat(rhs,P(L[i])*b[i])
        end
    end
    BasisExpansion(L[1].domain,Op\rhs)
end

function testconv(f::BasisExpansion)
    norm(f.c[end-4:end]) < tol
end


function \(L::Vector{ConcreteOperator},b::Vector)
    if !(typeof(N) <: Integer)
        n = 32
        sol = \(L,b,n)
        bool = testconv(sol)
        while !bool
            n *= 2
            sol = \(L,b,n)
            bool = testconv(sol)
        end
        return sol
    else
        return \(L,b,N)
    end
end

function \(L::Vector{AbstractOperator},b::Vector)
    Ops = [J*basis for J in L]
    Ops\b
end

function \(L::Vector{AbstractOperator},b::Vector,basis::Basis,N::Integer)
    Ops = [J*basis for J in L]
    \(Ops,b,N)
end

function \(L::Vector{AbstractOperator},b::Vector,N::Integer)
    Ops = [J*basis for J in L]
    \(Ops,b,N)
end

end # module OperatorApproximation