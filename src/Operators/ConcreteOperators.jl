abstract type MatrixOperator <: Operator end
#abstract type LazyOperator <: MatrixOperator end
abstract type BandedOperator <: MatrixOperator end  #likely don't need this distinction
abstract type SingleBandedOperator <: BandedOperator end

struct ConcreteOperator{D<:Basis,R<:Basis,T<:MatrixOperator} <: Operator
    domain::D
    range::R
    L::T
end

function *(Op::AbstractOperator,C::ConcreteOperator)  
    # Default.  Requires that the multiplication below
    # is defined...
    (Op*C.range)*C
end

function *(Op1::ConcreteOperator,Op2::ConcreteOperator)
    if Op1.domain != Op2.range
        @error "Domain-range mismatch."
    end
    ConcreteOperator(Op2.domain,Op1.range,Op1.L*Op2.L)
end

for op in (:+,:-)
    @eval begin 
        function $op(Op1::ConcreteOperator,Op2::ConcreteOperator) # need to check that the range and domain are compatible.
            #TODO: Check domain and range here
            if Op1.range != Op2.range
                @error "Range mismatch."
            elseif Op1.domain != Op2.domain
                @error "Domain mismatch."
            end
            ConcreteOperator(Op1.domain,Op1.range, $op(Op1.L,Op2.L))
        end
    end
end

function *(a::Number,Op::ConcreteOperator) 
    ConcreteOperator(Op.domain,Op.range, a*Op.L)
end

function Matrix(Op::ConcreteOperator{D,R,T},n,m) where {D,R,T}
    Matrix(Op.L,n,m)
end

function Matrix(Op::ConcreteOperator,m)
    Matrix(Op.L,m)
end

function *(Op::ConcreteOperator{D,R,T},f::BasisExpansion{S}) where {D,R,T,S}
    BasisExpansion(Op.range,Op.L*f.c)
end

function *(CC::Conversion,f::BasisExpansion)
    (CC*f.basis)*f
end

function rank(OP::ConcreteOperator)
    dim(OP.range)
end

struct SumOfMatrixOperators{T<:CoefficientDomain, S<:CoefficientDomain} <: MatrixOperator
    Ops::Vector{S} where S <: MatrixOperator
    c::Vector{Number} # be more specific?
end
SumOfMatrixOperators(Ops,c) = SumOfMatrixOperators{dom(Ops[1]),ran(Ops[1])}(Ops,c)

function *(a::Number,L::MatrixOperator)
    SumOfMatrixOperators([L],[a])
end

function *(a::Number,L::SumOfMatrixOperators)
    SumOfMatrixOperators(L.Ops,a*L.c)
end

function *(Op1::SumOfMatrixOperators,Op2::MatrixOperator)
    SumOfMatrixOperators([l*Op2 for l in Op1.Ops],Op1.c)
end

function getindex(Op::ConcreteOperator{D,R,T},n) where {D,R,T <: SumOfMatrixOperators}
    ConcreteOperator(Op.domain,Op.range,Op.L.Ops[n])
end

function Matrix(Op::ConcreteOperator{D,R,T},n,m) where {D,R,T <: SumOfMatrixOperators}
    A = Op.L.c[1]*Matrix(Op[1],n,m)
    for i = 2:length(Op.L.c)
        A += Op.L.c[i]*Matrix(Op[i],n,m)
    end
    A
end

function Matrix(Op::SumOfMatrixOperators,n,m)
    # TODO: check domain & range
    A = Op.c[1]*Matrix(Op.Ops[1],n,m)
    for i = 2:length(Op.c)
        A += Op.c[i]*Matrix(Op.Ops[i],n,m)
    end
    A
end

mutable struct SemiLazyBandedOperator{T<:CoefficientDomain, S<: CoefficientDomain} <: SingleBandedOperator# where T <: DiscreteDomain
    const nm::Integer
    const np::Integer
    const mat::Function
    A::Union{SparseMatrixCSC,Matrix}
end

struct BasicBandedOperator{T<:CoefficientDomain, S<: CoefficientDomain} <: SingleBandedOperator
    nm::Integer
    np::Integer
    A::Function
end

struct ZeroOperator{T<:CoefficientDomain, S<: CoefficientDomain} <: SingleBandedOperator
    nm::Integer
    np::Integer
    A::Function
end
ZeroOperator() = ZeroOperator{𝕏,𝕏}(0,0,x -> 0)

function *(Op::ZeroOperator,Op2::MatrixOperator)
    Op
end

function *(Op::MatrixOperator,Op2::ZeroOperator)
    Op2
end

struct ProductOfBandedOperators{T<:CoefficientDomain, S<:CoefficientDomain} <: BandedOperator
    V::Vector{S} where S <: SingleBandedOperator
end
ProductOfBandedOperators(V) = ProductOfBandedOperators{dom(V[end]),ran(V[1])}(V)

struct BlockMatrixOperator{T <: CoefficientDomain, S <: CoefficientDomain} <: MatrixOperator
    Ops::Matrix{MatrixOperator}
end
function BlockMatrixOperator(Ops)
    if typeof(Ops) <: Vector
        BlockMatrixOperator{𝕏,𝕏}(reshape(Ops,:,1))
    else
        BlockMatrixOperator{𝕏,𝕏}(Ops)
    end
end

# Repeat for +/- sign
# could be simplified using *
for op in (:+,:-)
    @eval begin
        function $op(L::SumOfMatrixOperators)
            SumOfLazyOperators(L.Ops,$op(L.c))
        end

        function $op(L::MatrixOperator)
            $op(1.0)*L
        end

        function $op(L1::SumOfMatrixOperators,L2::SumOfMatrixOperators)
            SumOfLazyOperators(vcat(L1.Ops,L2.Ops),vcat(L1.c,$op(L2.c)))
        end

        function $op(L1::SumOfMatrixOperators,L2::MatrixOperator)
            $op(L1,1.0*L2)
        end

        function $op(L1::MatrixOperator,L2::SumOfMatrixOperators)
            $op(1.0*L1,L2)
        end

        function $op(L1::MatrixOperator,L2::MatrixOperator)
            L1 + $op(1.0)*L2
        end

        function $op(L1::BlockMatrixOperator,L2::BlockMatrixOperator)
            BlockMatrixOperator($op.(L1.Ops,L2.Ops))
        end
    end
end

# Here S is not a DirectSum
function *(Op::ConcreteOperator{D,R,T},f::BasisExpansion{S}) where {D, R, T <: BlockMatrixOperator, S}
    BasisExpansion(Op.range,[Op.L.Ops[i,1]*f.c for i in 1:size(Op.L.Ops,1)])
end

# Here S is a DirectSum
function *(Op::ConcreteOperator{D,R,T},f::BasisExpansion{S}) where {D, R, T <: BlockMatrixOperator, S <: DirectSum}
    #return Op[1:end,1]
    out = Op[1:end,1]*f[1]
    for i = 2:size(Op,2)
        out += Op[1:end,i]*f[i]
    end
    out
end

size(L::BlockMatrixOperator) = size(L.Ops)
size(Op::ConcreteOperator{D,R,T}) where {D, R, T <: BlockMatrixOperator} = size(Op.L)

size(L::BlockMatrixOperator,j) = size(L.Ops,j)
size(Op::ConcreteOperator{D,R,T},j) where {D, R, T <: BlockMatrixOperator} = size(Op.L,j)

axes(L::BlockMatrixOperator) = axes(L.Ops)
axes(Op::ConcreteOperator{D,R,T}) where {D, R, T <: BlockMatrixOperator} = axes(Op.L)

axes(L::BlockMatrixOperator,j) = axes(L.Ops,j)
axes(Op::ConcreteOperator{D,R,T},j) where {D, R, T <: BlockMatrixOperator} = axes(Op.L,j)

function getindex(L::BlockMatrixOperator,i::Int64,j::Int64)
    L.Ops[i,j]
end
function getindex(L::BlockMatrixOperator,i::Int64,j::UnitRange{Int64})
    BlockMatrixOperator(L.Ops[i,j])
end
function getindex(L::BlockMatrixOperator,i::UnitRange{Int64},j::Int64)
    BlockMatrixOperator(L.Ops[i,j])
end
function getindex(L::BlockMatrixOperator,i::UnitRange{Int64},j::UnitRange{Int64})
    BlockMatrixOperator(L.Ops[i,j])
end

function getindex(Op::ConcreteOperator{D,R,T},i,j) where {D, R, T <: BlockMatrixOperator}
    ConcreteOperator(Op.domain[j],Op.range[i],Op.L[i,j])
end

# momtm = matrix of matrices to matrix
# will go recursive
function momtm(Ms::Union{Matrix{T},SparseMatrixCSC}) where T <: Union{Float64,ComplexF64}
    Ms
end

function momtm(Ms::Matrix{T}) where T <: Union{Matrix,SparseMatrixCSC,AbstractMatrix}
    A = hcat(momtm.(Ms[1,:])...)
    for i = 2:size(Ms,1)
        A = vcat(A,hcat(momtm.(Ms[i,:])...))
    end
    A
end

function Matrix(Op::BlockMatrixOperator,ns::Vector{Int64},ms::Vector{Int64})
    # Make robust to length of vectors?
    nmat = repeat(ns',length(ms))' |> Matrix
    mmat = repeat(ms',length(ns))
    momtm(Matrix.(Op.Ops,nmat,mmat))
end

# divides n in to k "equal-sized" bins
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

function divide_DOF(Op::ConcreteOperator{D,R,T},n::Integer,m::Integer) where {D <: Basis,R <: Basis,T <: BlockMatrixOperator}
    ns = divide_DOF(Op.range,n)
    ms = divide_DOF(Op.domain,m)
    ns, ms
end

function Matrix(Op::ConcreteOperator{D,R,T},n::Integer,m::Integer) where {D <: Basis,R <: Basis, T <: BlockMatrixOperator}
    ns, ms = divide_DOF(Op,n,m)
    Matrix(Op.L,ns,ms)
end

####
function ⊞(A1::MatrixOperator,A2::MatrixOperator)
    BlockMatrixOperator(reshape(A1,A2,1,:))
end

function ⊞(A1::BlockMatrixOperator,A2::MatrixOperator)
    BlockMatrixOperator(hcat(A1.Ops,A2))
end

function ⊞(A1::MatrixOperator,A2::BlockMatrixOperator)
    BlockMatrixOperator(hcat(A1,A2.Ops))
end

function ⊞(A1::BlockMatrixOperator,A2::BlockMatrixOperator)
    BlockMatrixOperator(hcat(A1.Ops,A2.Ops))
end
####
####
function ⊘(A1::MatrixOperator,A2::MatrixOperator)
    BlockMatrixOperator(reshape([A1 A2],:,1))
end

function ⊘(A1::BlockMatrixOperator,A2::MatrixOperator)
    BlockMatrixOperator(vcat(A1.Ops,A2))
end

function ⊘(A1::MatrixOperator,A2::BlockMatrixOperator)
    BlockMatrixOperator(vcat(A1,A2.Ops))
end

function ⊘(A1::BlockMatrixOperator,A2::BlockMatrixOperator)
    BlockMatrixOperator(vcat(A1.Ops,A2.Ops))
end
####
####

function ⊞(A1::ConcreteOperator,A2::ConcreteOperator)
    if A1.range == A2.range
        return ConcreteOperator(A1.domain ⊕ A2.domain, A1.range, A1.L ⊞ A2.L)
    else
        @error "Ranges are not compatible"
    end
end

function ⊘(A1::ConcreteOperator,A2::ConcreteOperator)
    if A1.domain == A2.domain
        return ConcreteOperator(A1.domain, A1.range ⊕ A2.range, A1.L ⊘ A2.L)
    else
        @error "Ranges are not compatible"
    end
end

dom(op::MatrixOperator) = collect(typeof(op).parameters)[1]
ran(op::MatrixOperator) = collect(typeof(op).parameters)[2]

for sop in (:BasicBandedOperator,:SemiLazyBandedOperator)
    @eval begin
        function *(A::$sop,B::$sop)
            ProductOfBandedOperators([A;B])
        end
    
        function *(A::$sop,B::ProductOfBandedOperators)
            ProductOfBandedOperators(vcat([A],B.V))
        end
    
        function *(B::ProductOfBandedOperators,A::$sop)
            ProductOfBandedOperators(vcat(B.V,[A]))
        end
    end
    for sop2 in (:BasicBandedOperator,:SemiLazyBandedOperator)
        if sop2 != sop
            @eval begin
                function *(A::$sop,B::$sop2)
                    ProductOfBandedOperators([A;B])
                end
            end
        end
    end
end
function *(B::ProductOfBandedOperators,A::ProductOfBandedOperators)
    ProductOfBandedOperators(vcat(B.V,A.V))
end

function *(Ops::SumOfMatrixOperators,Op::ProductOfBandedOperators)
    SumOfMatrixOperators([op*Op for op in Ops.Ops],Ops.c)
end

# include("SemiInfinite.jl")
# include("BiInfinite.jl")
include("LazyMatrix.jl")
include("Dense.jl")
include("GridValues/GridValues.jl")
include("Jacobi/Jacobi.jl")
include("Ultraspherical/Ultraspherical.jl")
include("Fourier/Fourier.jl")
include("Laurent/Laurent.jl")
include("Hardy/Hardy.jl")
include("DirectSum/DirectSum.jl")
include("Hermite/Hermite.jl")
include("Erf/Erf.jl")