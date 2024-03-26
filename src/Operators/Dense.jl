## Not really dense, per se, but no band structure

abstract type DenseOperator <: LazyOperator end
abstract type BasisEvaluationOperator <: DenseOperator end  # Always true for spectral methods
abstract type NaiveTransform <: DenseOperator end
abstract type FastTransform <: DenseOperator end

struct DenseTimesBanded{T <: CoefficientDomain, S <: CoefficientDomain} <: DenseOperator
    dense::TT where TT <: DenseOperator
    banded::SS where SS <: BandedOperator
end
DenseTimesBanded(dense,banded) = DenseTimesBanded{𝔼,dom(banded)}(dense,banded)

## This should resolve the ambiguity of the inner dimensions in the
## dense x dense multiplications.  But not now...
struct DenseTimesDense{T <: CoefficientDomain, S <: CoefficientDomain} <: DenseOperator
    denseL::DenseOperator
    denseR::DenseOperator
end
DenseTimesDense(denseL,denseR) = DenseTimesDense{𝔼,𝔼}(denseL,denseR)


function *(dense::DenseOperator,banded::BandedOperator)
    DenseTimesBanded(dense,banded)
end

function *(denseL::DenseOperator,denseR::DenseOperator)
    DenseTimesDense(denseL,denseR)
end

function Matrix(Op::DenseTimesBanded,n,m)
    nn = max(m + rowgrowth(Op.banded),0) # could be optimized
    B = Matrix(Op.banded,nn,m)
    Matrix(Op.dense,n,nn)*B
end

function Matrix(Op::DenseTimesDense,n,m)
    B = Matrix(Op.denseR,n,m)
    Matrix(Op.denseL,n,n)*B
end

function *(CC::Conversion,dom::Basis)
    if isconvertible(dom,CC.range)
        conversion(dom,CC.range) # convert from dom to CC.range
    else
        @error "Bases are not convertible."
    end
end

struct OPEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Function
    a::Function # Jacobi coefficients
    b::Function
end
OPEvaluationOperator(grid,a,b) = OPEvaluationOperator{ℕ₊,𝔼}(grid,a,b)

struct FourierEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Function
end
FourierEvaluationOperator(grid) = FourierEvaluationOperator{ℤ,𝔼}(grid)

struct FixedGridOPEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Vector
    a::Function # Jacobi coefficients
    b::Function
end
FixedGridOPEvaluationOperator(grid,a,b) = FixedGridOPEvaluationOperator{ℕ₊,𝔼}(grid,a,b)


struct FixedGridFourierEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Vector
end
FixedGridFourierEvaluationOperator(grid) = FixedGridFourierEvaluationOperator{ℤ,𝔼}(grid)

struct OPCauchyEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Function
    a::Function # Jacobi coefficients
    b::Function
    seed::Function
end
OPCauchyEvaluationOperator(grid,a,b,seed) = OPCauchyEvaluationOperator{ℕ₊,𝔼}(grid,a,b,seed)

mutable struct OPEigenTransform{T <: CoefficientDomain, S <: CoefficientDomain} <: NaiveTransform
    const a::Function # Jacobi coefficients
    const b::Function
    A::Matrix # saves transform matrix
    function OPEigenTransform{𝔼,ℕ₊}(a,b)
        return new(a,b,hcat(1.0))
    end
end
OPEigenTransform(a,b) = OPEigenTransform{𝔼,ℕ₊}(a,b)

struct DiscreteFourierTransform{T <: CoefficientDomain, S <: CoefficientDomain} <: FastTransform 
    T::Function
    function DiscreteFourierTransform{𝔼,ℤ}()
        return new(mdft)
    end
end
DiscreteFourierTransform() = DiscreteFourierTransform{𝔼,ℤ}()

struct GridMultiplication{T <: CoefficientDomain, S <: CoefficientDomain} <: DenseOperator # even though it is sparse...
    # it is simpler to treat grid multiplication as dense
    f::Function
    grid::Function
end
GridMultiplication(f,grid) = GridMultiplication{𝔼,𝔼}(f,grid)

struct FixedGridMultiplication{T <: CoefficientDomain, S <: CoefficientDomain} <: DenseOperator
    fvals::Vector
end
FixedGridMultiplication(f,grid) = FixedGridMultiplication{𝔼,𝔼}(f,grid)

function Matrix(Op::OPEvaluationOperator,n,m)
    poly(Op.a,Op.b,m,Op.grid(n)) 
end

function horner_mat(x,m)
    A = zeros(ComplexF64,length(x),m)
    mm = convert(Int64,floor( m/2 ))
    A[:,1] = exp.(-1im*pi*mm*x)
    ex1 = exp.(1im*pi*x)
    for i = 2:m
        A[:,i]  =  copy(A[:,i-1]).*ex1
    end
    return A
end

function Matrix(Op::FourierEvaluationOperator,n,m)
    hornermat(Op.grid(n),m)
end

function Matrix(Op::FixedGridOPEvaluationOperator,n,m)
    if n <= length(Op.grid)
        return poly(Op.a,Op.b,m,Op.grid[1:n])
    else
        @warn "Asked for more rows than grid points.  Returning maximum number of rows."
        return poly(Op.a,Op.b,m,Op.grid)
    end
end

function Matrix(Op::FixedGridFourierEvaluationOperator,n,m)
    if n <= length(Op.grid)
        return hornermat(Op.grid[1:n],m)
    else
        @warn "Asked for more rows than grid points.  Returning maximum number of rows."
        return hornermat(Op.grid,m)
    end
end

function Matrix(Op::FixedGridOPEvaluationOperator,m)  # only one dim for Functional
    return poly(Op.a,Op.b,m,Op.grid)
end

function Matrix(Op::OPCauchyEvaluationOperator,n,m)
    cauchy(Op.a,Op.b,seed,m,Op.grid(n)) 
end

function Matrix(Op::OPEigenTransform,n)
    if size(Op.A)[1] == n
        return Op.A
    end
    Op.A = Interp_transform(Op.a,Op.b,n-1)[2]
    return Op.A
end

function Matrix(Op::OPEigenTransform,n,m)
    if size(Op.A)[1] != m
        Op.A = Interp_transform(Op.a,Op.b,m-1)[2]
    end
    if n == m
        return Op.A
    elseif n < m
        return Op.A[1:n,1:m]
    else
        return vcat(Op.A,zeros(n-m,m))
    end
end

# TODO: Use Clenshaw
function *(Op::DenseOperator,v::Vector)
    Matrix(Op,length(v))*v
end

function *(Op::FastTransform,v::Vector)
    Op.T(v)
end

function Matrix(Op::DiscreteFourierTransform,n,m)
    Op.T(Matrix(I,n,m)) # Not the right way to do this...
end
