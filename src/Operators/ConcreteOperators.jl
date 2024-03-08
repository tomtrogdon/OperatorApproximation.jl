abstract type LazyOperator <: Operator end
abstract type ConcreteOperator <: Operator end

struct ConcreteLazyOperator{D<:Basis,R<:Basis,T<:LazyOperator} <: ConcreteOperator
    domain::D
    range::R
    L::T
end

function *(Op::AbstractOperator,C::ConcreteLazyOperator)  
    # Default.  Requires that the multiplication below
    # is defined...
    (Op*C.range)*C
end

function *(Op1::ConcreteLazyOperator,Op2::ConcreteLazyOperator)
    if Op1.domain != Op2.range
        @error "Domain-range mismatch."
    end
    ConcreteLazyOperator(Op2.domain,Op1.range,Op1.L*Op2.L)
end

for op in (:+,:-)
    @eval begin 
        function (Op1::ConcreteLazyOperator,Op2::ConcreteLazyOperator) # need to check that the range and domain are compatible.
            #TODO: Check domain and range here
            if Op1.range != Op2.range
                @error "Range mismatch."
            elseif Op1.domain != Op2.domain
                @error "Domain mismatch."
            end
            ConcreteLazyOperator(Op1.domain,Op1.range, $op(Op1.L,Op2.L))
        end
    end
end

function *(a::Number,Op::ConcreteLazyOperator) 
    ConcreteLazyOperator(Op.domain,Op.range, a*Op.L)
end

function Matrix(Op::ConcreteLazyOperator{D,R,T},n,m) where {D,R,T}
    Matrix(Op.L,n,m)
end

function Matrix(Op::ConcreteLazyOperator,m)
    Matrix(Op.L,m)
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
SI = NN()  ## Semi-Infinite
BI = ZZ()  ## Bi-Infinite

abstract type BandedOperator <: LazyOperator end
abstract type SingleBandedOperator <: BandedOperator end


## Note: Needs DiscreteDomain field if products of sums are to be used
struct SumOfLazyOperators <: LazyOperator
    Ops::Vector{S} where S <: LazyOperator
    c::Vector{Number} # be more specific?
end

function *(a::Number,L::LazyOperator)
    SumOfLazyOperators([L],[a])
end

function *(a::Number,L::SumOfLazyOperators)
    SumOfLazyOperators(L.Ops,a*L.c)
end

# Repeat for +/- sign
# could be simplified using *
for op in (:+,:-)
    @eval begin
        function $op(L::SumOfLazyOperators)
            SumOfLazyOperators(L.Ops,$op(L.c))
        end

        function $op(L::LazyOperator)
            $op(1.0)*L
        end

        function $op(L1::SumOfLazyOperators,L2::SumOfLazyOperators)
            SumOfLazyOperators(vcat(L1.Ops,L2.Ops),vcat(L1.c,$op(L2.c)))
        end

        function $op(L1::SumOfLazyOperators,L2::LazyOperator)
            $op(L1,1.0*L2)
        end

        function $op(L1::LazyOperator,L2::SumOfLazyOperators)
            $op(1.0*L1,L2)
        end

        function $op(L1::LazyOperator,L2::LazyOperator)
            L1 + $op(1.0)*L2
        end
    end
end

function Matrix(Op::SumOfLazyOperators,n,m)
    # TODO: check domain & range
    A = Op.c[1]*Matrix(Op.Ops[1],n,m)
    for i = 2:length(Op.c)
        A += Op.c[i]*Matrix(Op.Ops[i],n,m)
    end
    A
end

mutable struct SemiLazyBandedOperator{T<:DiscreteDomain} <: SingleBandedOperator# where T <: DiscreteDomain
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

struct BlockLazyOperator <: LazyOperator
    Ops::Matrix{LazyOperator}
end

# momtm = matrix of matrices to matrix
# will go recursive
function momtm(Ms::Union{Matrix{T},SparseMatrixCSC}) where T <: Union{Float64,ComplexF64}
    Ms
end

function momtm(Ms::Matrix{T}) where T <: Union{Matrix,SparseMatrixCSC}
    A = hcat(momtm.(Ms[1,:])...)
    for i = 2:size(Ms,1)
        A = vcat(A,hcat(momtm.(Ms[i,:])...))
    end
    A
end

function Matrix(Op::BlockLazyOperator,ns::Vector{Int64},ms::Vector{Int64})
    # Make robust to length of vectors?
    nmat = repeat(ns',length(ms))' |> Matrix
    mmat = repeat(ms',length(ns))
    momtm(Matrix.(Op.Ops,nmat,mmat))
end

function binit(n,k)
    ((n+k-1:-1:n) .÷ k) |> Vector
end

function divide_DOF(b::Basis,n::Integer)
    N = n
    ranges = bases(b)
    ns = zeros(Int64,length(ranges))
    for i = 1:length(ranges)
        if dim(ranges[i]) < Inf
            ns[i] = dim(ranges[i])
        end
    end
    N -= sum(ns)
    if N < 1
        @error "n, m not large enough to capture functionals"
        return
    end
    ninds = ns .== 0
    ns[ninds] = binit(N, sum(ninds))
    ns
end

function divide_DOF(Op::ConcreteLazyOperator{D,R,T},n::Integer,m::Integer) where {D <: Basis,R <: Basis,T <: BlockLazyOperator}
    ns = divide_DOF(Op.range,n)
    ms = divide_DOF(Op.domain,m)
    ns, ms
end

function Matrix(Op::ConcreteLazyOperator{D,R,T},n::Integer,m::Integer) where {D <: Basis,R <: Basis,T <: BlockLazyOperator}
    ns, ms = divide_DOF(Op,n,m)
    Matrix(Op.L,ns,ms)
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
        for sop2 in (:BasicBandedOperator,:SemiLazyBandedOperator)
            if sop2 != sop
                @eval begin
                    function *(A::$sop{T},B::$sop2{S}) where {S <: $op, T<: $op}
                        MultipliedBandedOperator(A.DD,[A;B])
                    end
                end
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