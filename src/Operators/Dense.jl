## Not really dense, per se, but no band structure

abstract type DenseOperator <: MatrixOperator end  #sometimes P
abstract type BasisEvaluationOperator <: DenseOperator end  # Always true for spectral methods
abstract type NaiveTransform <: DenseOperator end
abstract type FastTransform <: DenseOperator end

struct FixedMatrix{T<:CoefficientDomain, S <: CoefficientDomain} <: DenseOperator# where T <: DiscreteDomain
    A::Union{SparseMatrixCSC,Matrix}
end
FixedMatrix(A) = FixedMatrix{𝕏,𝕏}(A) # what should default be?

function Matrix(Op::FixedMatrix,n,m)
    nn, mm = size(Op.A)
    if n > nn
        @warn "MatrixOperator: n is too large"
    end
    if m > mm
        @warn "MatrixOperator: m is too large"
    end
    Op.A[1:min(n,nn),1:min(m,mm)]
end

struct DenseTimesBanded{T <: CoefficientDomain, S <: CoefficientDomain} <: DenseOperator
    dense::TT where TT <: DenseOperator
    banded::SS where SS <: BandedOperator
end
DenseTimesBanded(dense,banded) = DenseTimesBanded{ran(dense),dom(banded)}(dense,banded)

## This should resolve the ambiguity of the inner dimensions in the
## dense x dense multiplications.  But not now...
struct DenseTimesDense{T <: CoefficientDomain, S <: CoefficientDomain} <: DenseOperator
    denseL::DenseOperator
    denseR::DenseOperator
end
DenseTimesDense(denseL,denseR) = DenseTimesDense{ran(denseL),dom(denseR)}(denseL,denseR)

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

function *(CC::CoefConversion,dom::Basis)
    if cfd(dom) == cfd(CC.range)
        ConcreteOperator(dom,CC.range,BasicBandedOperator{cfd(dom),cfd(dom)}(0,0,(i,j) ->  i == j ? 1.0 : 0.0))
         # convert from dom to CC.range
    else
        @error "Bases are not coef-convertible."
    end
end

struct OPEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Union{Function,Vector}
    a::Function # Jacobi coefficients
    b::Function
end
OPEvaluationOperator(grid,a,b) = OPEvaluationOperator{ℕ₊,𝔼}(grid,a,b)

struct WeightedOPEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Union{Function,Vector}
    a::Function # Jacobi coefficients
    b::Function
    W::Function
end
WeightedOPEvaluationOperator(grid,a,b) = WeightedOPEvaluationOperator{ℕ₊,𝔼}(grid,a,b)

struct FourierEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Union{Function,Vector}
end
FourierEvaluationOperator(grid) = FourierEvaluationOperator{ℤ,𝔼}(grid)

struct RationalEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Union{Function,Vector}
end
RationalEvaluationOperator(grid) = RationalEvaluationOperator{ℤ,𝔼}(grid)

struct OPCauchyEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Union{Function,Vector}
    a::Function # Jacobi coefficients
    b::Function
    seed::Function
end
OPCauchyEvaluationOperator(grid,a,b,seed) = OPCauchyEvaluationOperator{ℕ₊,𝔼}(grid,a,b,seed)

struct PoleResCauchyEvaluationOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: BasisEvaluationOperator
    grid::Union{Function,Vector}
    ps::Vector # poles
end
PoleResCauchyEvaluationOperator(grid,ps) = PoleResCauchyEvaluationOperator{ℕ₊,𝔼}(grid,ps)

mutable struct OPEigenTransform{T <: CoefficientDomain, S <: CoefficientDomain} <: NaiveTransform
    const a::Function # Jacobi coefficients
    const b::Function
    A::Matrix # saves transform matrix
    function OPEigenTransform{𝔼,ℕ₊}(a,b)
        return new(a,b,hcat(1.0))
    end
end
OPEigenTransform(a,b) = OPEigenTransform{𝔼,ℕ₊}(a,b)

mutable struct OPWeightedEigenTransform{T <: CoefficientDomain, S <: CoefficientDomain} <: NaiveTransform
    const a::Function # Jacobi coefficients
    const b::Function
    A::Matrix # saves transform matrix
    W::Function
    function OPWeightedEigenTransform{𝔼,ℕ₊}(a,b,W)
        return new(a,b,hcat(1.0),W)
    end
end
OPWeightedEigenTransform(a,b,W) = OPWeightedEigenTransform{𝔼,ℕ₊}(a,b,W)

struct DiscreteFourierTransform{T <: CoefficientDomain, S <: CoefficientDomain} <: FastTransform 
    T::Function
    function DiscreteFourierTransform{𝔼,ℤ}()
        return new(mdft)
    end
end
DiscreteFourierTransform() = DiscreteFourierTransform{𝔼,ℤ}()

function *(D::DiscreteFourierTransform,f::Vector)
    D.T(f)
end

struct GridMultiplication{T <: CoefficientDomain, S <: CoefficientDomain} <: DenseOperator # even though it is sparse...
    # it is simpler to treat grid multiplication as dense
    f::Union{Function,Vector}
    grid::Union{Function,Vector}
end
GridMultiplication(f,grid) = GridMultiplication{𝔼,𝔼}(f,grid)

function Matrix(Op::OPEvaluationOperator,n,m)
    if typeof(Op.grid) <: Function
        return poly(Op.a,Op.b,m,Op.grid(n)) 
    end
    if n <= length(Op.grid)
        return poly(Op.a,Op.b,m,Op.grid[1:n])
    else
        @warn "Asked for more rows than grid points.  Returning maximum number of rows."
        return poly(Op.a,Op.b,m,Op.grid)
    end
end

function Matrix(Op::WeightedOPEvaluationOperator,n,m)
    if typeof(Op.grid) <: Function
        return Diagonal(Op.W.(Op.grid(n)))*poly(Op.a,Op.b,m,Op.grid(n)) 
    end
    if n <= length(Op.grid)
        return Diagonal(Op.W.(Op.grid[1:n]))*poly(Op.a,Op.b,m,Op.grid[1:n])
    else
        @warn "Asked for more rows than grid points.  Returning maximum number of rows."
        return Diagonal(Op.W.(Op.grid))*poly(Op.a,Op.b,m,Op.grid)
    end
end

function horner_mat(x,m)
    A = zeros(ComplexF64,length(x),m)
    mm = convert(Int64,floor( m/2 ))
    A[:,1] = exp.(-1im*pi*mm*x)
    ex1 = exp.(1im*pi*x)
    for i = 2:m
        A[:,i]  .=  copy(A[:,i-1]).*ex1
    end
    return A
end

function Matrix(Op::FourierEvaluationOperator,n,m)
    if typeof(Op.grid) <: Function
        return hornermat(Op.grid(n),m)
    end
    if n <= length(Op.grid)
        return hornermat(Op.grid[1:n],m)
    else
        @warn "Asked for more rows than grid points.  Returning maximum number of rows."
        return hornermat(Op.grid,m)
    end
end

function horner_mat_rat(x,m)
    #need to ensure that the x that's being passed in is equivalent to:
    #     x = P.basis.GD.D.imap(k)
    #     x = ((x.-1im)./(x.+1im))
    A = zeros(ComplexF64,length(x),m)
    mm = convert(Int64,floor( m/2 ))
    A[:,1] = x^(-mm)
    ex1 = x^(-mm)
    for i = 2:m
        A[:,i]  .=  copy(A[:,i-1]).*ex1
    end
    return A
end

function Matrix(Op::RationalEvaluationOperator,n,m,α)
    if typeof(Op.grid) <: Function
        k = Op.grid(n) #is this a row vector or a column vector? Because I may need to transpose this in exp term
        x = ((k.-1im)./(k.+1im))
        return horner_mat_rat(x,m).*exp(1im*k*α)
    end
    if n <= length(Op.grid)
        k = Op.grid[1:n]
        x = ((k.-1im)./(k.+1im))
        return horner_mat_rat(x,m).*exp(1im*k*α)
    else
        @warn "Asked for more rows than grid points.  Returning maximum number of rows."
        k = Op.grid
        x = ((k.-1im)./(k.+1im))
        return horner_mat_rat(x,m).*exp(1im*k*α)
    end
end

function Matrix(Op::OPEvaluationOperator,m)  # only one dim for Functional
    if typeof(Op.grid) <: Function
        @error "Need another dimension"
        return
    end
    return poly(Op.a,Op.b,m,Op.grid)
end

function Matrix(Op::OPCauchyEvaluationOperator,n,m)
    if typeof(Op.grid) <: Function
        return cauchy(Op.a,Op.b,Op.seed,m-1,Op.grid(n))*2
    else
        return return cauchy(Op.a,Op.b,Op.seed,m-1,Op.grid)*2
    end
end

function Matrix(Op::PoleResCauchyEvaluationOperator,n,m)
    if m != length(Op.ps)
        @warn "PoleResCauchyEvaluationOperator: Incorrect residue dim"
    end
    if typeof(Op.grid) <: Function
        return poleres_cauchy(Op.ps,Op.grid(n))
    else
        if n > length(Op.grid)
            @warn "PoleResCauchyEvaluationOperator: Incorrect grid dim"
            return poleres_cauchy(Op.ps,Op.grid)
        end
        return poleres_cauchy(Op.ps,Op.grid[1:n])
    end
end

function Matrix(Op::OPEigenTransform,n)
    if size(Op.A)[1] == n
        return Op.A
    end
    Op.A = Interp_transform(Op.a,Op.b,n-1)[2]
    return Op.A
end

function Matrix(Op::OPWeightedEigenTransform,n)
    if size(Op.A)[1] == n
        return Op.A
    end
    λ, O = Interp_transform(Op.a,Op.b,n-1)
    Op.A = O*Diagonal(Op.W(λ))
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

function Matrix(Op::OPWeightedEigenTransform,n,m)
    if size(Op.A)[1] != m
        λ, O = Interp_transform(Op.a,Op.b,n-1)
        Op.A = O*Diagonal(Op.W(λ))
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
