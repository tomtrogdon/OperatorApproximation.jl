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

function *(Op::SumOfAbstractOperators,C::ConcreteLazyOperator)
    if Op.domain != C.range
        @error "Domain-range mismatch."
    end
    ops = [op*C for op in Op.Ops]
    SumOfConcreteOperators(ops[1].domain,ops[1].range,ops,Op.c)
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

function *(CC::Conversion,f::BasisExpansion)
    (CC*f.basis)*f
end

function rank(OP::ConcreteOperator)
    dim(OP.range)
end


abstract type DiscreteDomain end
struct ZZ <: DiscreteDomain end
struct NN <: DiscreteDomain end
SI = NN()
BI = ZZ()

abstract type BandedOperator <: LazyOperator end
abstract type SingleBandedOperator <: BandedOperator end

mutable struct SemiLazyBandedOperator{T<:DiscreteDomain} <: SingleBandedOperator where T <: DiscreteDomain
    const DD::T
    const nm::Integer
    const np::Integer
    const mat::Function
    A::SparseMatrixCSC
end

struct BasicBandedOperator{T<:DiscreteDomain} <: SingleBandedOperator
    DD::T
    nm::Integer
    np::Integer
    A::Function
end

struct MultipliedBandedOperator{T<:DiscreteDomain} <: BandedOperator
    DD::T
    V::Vector{S} where S <: SingleBandedOperator
end

for op in (:ZZ,:NN)
    for sop in (:BasicBandedOperator,:SemiLazyBandedOperator)
        @eval begin
            function *(A::$sop{T},B::$sop{S}) where {S <: $op, T<: $op}
                MultipliedBandedOperator(A.DD,[A;B])
            end
    
            function *(A::$sop{T},B::MultipliedBandedOperator{S}) where {S <: $op, T<: $op}
                MultipliedBandedOperator(A.DD,vcat([A],B.V))
            end
    
            function *(B::MultipliedBandedOperator{T},A::$sop{S}) where {S <: $op, T<: $op}
                MultipliedBandedOperator(A.DD,vcat(B.V,[A]))
            end
        end
    end
    @eval begin
        function *(B::MultipliedBandedOperator{T},A::MultipliedBandedOperator{S}) where {S <: $op, T<: $op}
            MultipliedBandedOperator(A.DD,vcat(B.V,A.V))
        end
    end
end

include("SemiInfinite.jl")
include("BiInfinite.jl")
include("Dense.jl")

include("GridValues/GridValues.jl")
include("Jacobi/Jacobi.jl")
include("Ultraspherical/Ultraspherical.jl")
include("Fourier/Fourier.jl")