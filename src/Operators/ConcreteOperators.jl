abstract type LazyOperator <: Operator end
abstract type ConcreteOperator <: Operator end

struct ConcreteLazyOperator{D<:Basis,R<:Basis} <: ConcreteOperator
    domain::D
    range::R
    L::LazyOperator
end

struct SumOfConcreteOperators{D<:Basis,R<:Basis,T} <: ConcreteOperator where T
    domain::D
    range::R
    Ops::Vector{T}
    c::Vector
end

function *(Op::AbstractOperator,C::ConcreteLazyOperator)  
    # Default.  Requires that the multiplication below
    # is defined...
    # println("default")
    (Op*C.range)*C
end

function *(Op::SumOfAbstractOperators,C::ConcreteLazyOperator)  #TODO: Enforce domain/range
    if Op.domain != C.range
        @error "Domain-range mismatch."
    end
    ops = [op*C for op in Op.Ops]
    SumOfConcreteOperators(ops[1].domain,ops[1].range,ops,Op.c)
end

abstract type BandedOperator <: LazyOperator end
abstract type DenseOperator <: LazyOperator end
abstract type BasisEvaluationOperator <: DenseOperator end  # Always true for spectral methods

struct SingleBandedOperator <: BandedOperator
    nm::Integer
    np::Integer
    A::Function
end

function *(Op1::ConcreteLazyOperator,Op2::ConcreteLazyOperator)
    if Op1.domain != Op2.range
        @error "Domain-range mismatch."
    end
    ConcreteLazyOperator(Op2.domain,Op1.range,Op1.L*Op2.L)
end

function +(Op1::ConcreteLazyOperator,Op2::ConcreteLazyOperator) # need to check that the range and domain are compatible.
    #TODO: Check domain and range here
    if Op1.range != Op2.range
        @error "Range mismatch."
    elseif Op1.domain != Op2.domain
        @error "Domain mismatch."
    end
    SumOfConcreteOperators(Op1.domain,Op1.range,[Op1;Op2],[1;1])
end

 function Matrix(Op::SingleBandedOperator,n,m)
    A = spzeros(n,m)
    if Op.nm > n
        nm = n
    else
        nm = Op.nm
    end
    
    if Op.np > m
        np = m
    else
        np = Op.np
    end
    
    for i = -nm:-1
        if -i <= n-1
            A += spdiagm(n,m, (i => [Op.A(1 -i + j,1 + j) for j in 0:min(n+i-1,m-1) ]))
        end
    end
    
    for i = 1:np
        if i <= m-1
            A += spdiagm(n,m, (i => [Op.A(1 + j,1 + j + i) for j in 0:min(m-i-1,n-1)]))
        end
    end
    
    if Op.nm >= 0 && Op.np >= 0
        A += spdiagm(n,m, [Op.A(i,i) for i in 1:min(n,m)] )
    end
    A
end

struct MultipliedBandedOperator <: BandedOperator
    V::Vector{SingleBandedOperator}
end

function *(A::SingleBandedOperator,B::SingleBandedOperator)
    MultipliedBandedOperator([A;B])
end

function *(A::SingleBandedOperator,B::MultipliedBandedOperator)
    MultipliedBandedOperator(vcat([A],B.V))
end
    
function *(B::MultipliedBandedOperator,A::SingleBandedOperator)
    MultipliedBandedOperator(vcat(B.V,[A]))
end
        
function *(B::MultipliedBandedOperator,A::MultipliedBandedOperator)
    MultipliedBandedOperator(vcat(B.V,A.V))
end

function Matrix(Op::MultipliedBandedOperator,n,m)
    cols = m
    rows = max(cols+Op.V[end].nm,1)
    A = Matrix(Op.V[end],rows,cols)
    for j = length(Op.V)-1:-1:2
        cols = rows
        rows = max(cols + Op.V[j].nm,1)
        A = Matrix(Op.V[j],rows,cols)*A
    end
    cols = rows
    rows = n
    A = Matrix(Op.V[1],n,cols)*A
    A
end

function *(Op::SingleBandedOperator,c::Vector)
    n = length(c)
    m = max(Op.nm + n,1)
    Matrix(Op,m,n)*c
end

function *(Op::MultipliedBandedOperator,c::Vector)
    cols = length(c)
    rows = max(cols+Op.V[end].nm,1)
    v = Matrix(Op.V[end],rows,cols)*c
    for j = length(Op.V)-1:-1:1
        cols = rows
        rows = max(cols + Op.V[j].nm,1)
        v = Matrix(Op.V[j],rows,cols)*v
    end
    v
end

function rowgrowth(Op::SingleBandedOperator)
    Op.nm
end

function rowgrowth(Op::MultipliedBandedOperator)
    cols = 0
    rows = cols+Op.V[end].nm
    for j = length(Op.V)-1:-1:1
        cols = rows
        rows = cols + Op.V[j].nm
    end
    rows
end

function Matrix(Op::ConcreteLazyOperator,n,m)
    Matrix(Op.L,n,m)
end

function Matrix(Op::ConcreteLazyOperator,m)
    Matrix(Op.L,m)
end

function Matrix(Op::SumOfConcreteOperators,n,m)
    # TODO: check domain & range
    A = Op.c[1]*Matrix(Op.Ops[1],n,m)
    for i = 2:length(Op.c)
        A += Op.c[i]*Matrix(Op.Ops[i],n,m)
    end
    A
end

function *(Op::ConcreteLazyOperator,f::BasisExpansion)
    BasisExpansion(Op.range,Op.L*f.c)
end


struct DenseTimesBanded <: DenseOperator
    banded::BandedOperator
    dense::DenseOperator
end

## This should resolve the ambiguity of the inner dimensions in the
## dense x dense multiplications.  But not now...
struct DenseTimesDense <: DenseOperator
    denseL::DenseOperator
    denseR::DenseOperator
end

function *(dense::DenseOperator,banded::BandedOperator)
    DenseTimesBanded(banded,dense)
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

function *(CC::Conversion,f::BasisExpansion)
    (CC*f.basis)*f
end

struct OPEvaluationOperator <: BasisEvaluationOperator  ## Add CollocatedOperator?
    grid::Function
    a::Function # Jacobi coefficients
    b::Function
end

struct FixedGridOPEvaluationOperator <: BasisEvaluationOperator  ## Add CollocatedOperator?
    grid::Vector
    a::Function # Jacobi coefficients
    b::Function
end

function Matrix(Op::OPEvaluationOperator,n,m)
    poly(Op.a,Op.b,m,Op.grid(n)) 
end

function Matrix(Op::FixedGridOPEvaluationOperator,n,m)
    if n <= length(Op.grid)
        return poly(Op.a,Op.b,m,Op.grid[1:n])
    else
        @warn "Asked for more rows than grid points.  Returning maximum number of rows."
        return poly(Op.a,Op.b,m,Op.grid)
    end
end

function Matrix(Op::FixedGridOPEvaluationOperator,m)  # only one dim for Functional
    return poly(Op.a,Op.b,m,Op.grid)
end

function rank(OP::ConcreteOperator)
    dim(OP.range)
end

## Need to figure out operator Multiplication here.

include("GridValues/GridValues.jl")
include("Jacobi/Jacobi.jl")
include("Ultraspherical/Ultraspherical.jl")