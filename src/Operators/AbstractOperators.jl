#### SETTING UP A NEW OPERATOR ####
# (1) Set up an abstract operator type.   If the range will vary based on the setting,
#     have that be the only field.
# (2) Within the directory for a given basis, create a file for this operator.
#     Within this file describe how the product of the abstract operator should behave
#     on the basis.  The result should be a ConcreteOperator, with domain and range identified.
#     It could be the case that the product of this ConcreteOperator and another is not defined
#     or is ambiguious.  In this case, new ConcreteOperator structs may have be written.
abstract type Operator end

abstract type AbstractOperator <: Operator end

function *(Op::AbstractOperator,V::Vector)
    [Op*v for v in V]
end

function *(Op::AbstractOperator,V::Matrix)
    [Op*v for v in V]
end

# Note that we don't allow blocks of blocks (but maybe should...)
struct BlockAbstractOperator{T} <: AbstractOperator where T <: AbstractOperator
    Ops::Matrix{T}
end
function BlockAbstractOperator(Op::AbstractOperator,n,m)
    if n == m && n == 1
        return Op
    end
    Ops = fill(Op,n,m)
    BlockAbstractOperator(Ops)
end
function matrix2BlockOperator(M::Matrix{T}) where T <: Operator
    n,m = size(M)
    row = M[1,1]
    for j = 2:m
        row = row ⊞ M[1,j]
    end
    op = row
    for i = 2:n
        row = M[i,1]
        for j = 2:m
            row = row ⊞ M[i,j]
        end
        op = op ⊘ row
    end
    op
end

function diagm(V::Vector{T}) where T <: AbstractOperator
    if length(V) == 1
        return V[1]
    end
    c = convert(Vector{Int64},[])
    r = convert(Vector{Int64},[])
    for v in V
        n,m = size(v)
        push!(c,m)
        push!(r,n)
    end
    row = [AbstractZeroOperator(r[1],c[j]) for j in 1:length(V)]
    row = convert(Vector{AbstractOperator},row)
    row[1] = V[1]
    out = ⊞(row...)
    for i = 2:length(V)
        row = [AbstractZeroOperator(r[i],c[j]) for j in 1:length(V)]
        row = convert(Vector{AbstractOperator},row)
        row[i] = V[i]
        out = out ⊘ (⊞(row...))
    end
    out
end

struct BlockDiagonalAbstractOperator{T} <: AbstractOperator where T <: AbstractOperator
    Ops::Vector{T}
end
diagm(Op::BlockDiagonalAbstractOperator) = diagm(Op.Ops)

size(B::AbstractOperator) = (1,1)
size(B::BlockAbstractOperator) = size(B.Ops)
size(B::BlockAbstractOperator,i) = size(B.Ops,i)
size(B::BlockDiagonalAbstractOperator) = size(B.Ops)
axes(B::BlockAbstractOperator,j) = axes(B.Ops,j)

function getindex(L::BlockAbstractOperator,i::Int64,j::Int64)
    L.Ops[i,j]
end
function getindex(L::BlockAbstractOperator,i::Int64,j::UnitRange{Int64})
    BlockAbstractOperator(reshape(L.Ops[i,j],1,:))
end
function getindex(L::BlockAbstractOperator,i::UnitRange{Int64},j::Int64)
    BlockAbstractOperator(reshape(L.Ops[i,j],:,1))
end
function getindex(L::BlockAbstractOperator,i::UnitRange{Int64},j::UnitRange{Int64})
    BlockAbstractOperator(L.Ops[i,j])
end

####
function ⊕(A1::AbstractOperator,A2::AbstractOperator)
    BlockDiagonalAbstractOperator([A1, A2])
end

function ⊕(A1::BlockDiagonalAbstractOperator,A2::AbstractOperator)
    BlockDiagonalAbstractOperator(vcat(A1.Ops, [A2]))
end

function ⊕(A1::AbstractOperator,A2::BlockDiagonalAbstractOperator)
    BlockDiagonalAbstractOperator(vcat([A1], A2.Ops))
end

function ⊕(A1::BlockDiagonalAbstractOperator,A2::BlockDiagonalAbstractOperator)
    BlockDiagonalAbstractOperator(vcat(A1.Ops, A2.Ops))
end
####

####
function ⊞(A1::AbstractOperator,A2::AbstractOperator)
    BlockAbstractOperator([A1 A2])
end

function ⊞(A1::BlockAbstractOperator,A2::AbstractOperator)
    BlockAbstractOperator([A1.Ops A2])
end

function ⊞(A1::AbstractOperator,A2::BlockAbstractOperator)
    BlockAbstractOperator([A1 A2.Ops])
end

function ⊞(A1::BlockAbstractOperator,A2::BlockAbstractOperator)
    BlockAbstractOperator([A1.Ops A2.Ops])
end

function ⊞(Ops...)
    if length(Ops) == 1
        return Ops
    end
    op = Ops[1]
    for i = 2:length(Ops)
        op = op ⊞ Ops[i]
    end
    op
end

####
####
function ⊘(A1::AbstractOperator,A2::AbstractOperator)
    BlockAbstractOperator(reshape([A1 A2],:,1))
end

function ⊘(A1::BlockAbstractOperator,A2::AbstractOperator)
    BlockAbstractOperator([A1.Ops; A2])
end

function ⊘(A1::AbstractOperator,A2::BlockAbstractOperator)
    BlockAbstractOperator([A1; A2.Ops])
end

function ⊘(A1::BlockAbstractOperator,A2::BlockAbstractOperator)
    BlockAbstractOperator([A1.Ops; A2.Ops])
end
####
function *(D::BlockDiagonalAbstractOperator,B::BlockAbstractOperator)
    if length(D.Ops) != size(B,1)
        @error "Dimensions are incorrect in BlockDiagonalAbstractOperator * BlockAbstractOperator"
    end
    Ops = reshape([D.Ops[1]*b for b in B.Ops[1,1:end]],1,:)
    for i in 2:length(D.Ops)
        Ops = vcat(Ops,reshape([D.Ops[i]*b for b in B.Ops[i,1:end]],1,:))
    end
    BlockAbstractOperator(Ops)
end

function *(D::BlockAbstractOperator,B::BlockAbstractOperator)
    BlockAbstractOperator(convert(Matrix{OperatorApproximation.AbstractOperator},D.Ops*B.Ops))
end

function *(D::AbstractOperator,B::BlockAbstractOperator)
    if size(B.Ops,1) != 1
        @error "wrong block size"
    end
    BlockAbstractOperator([D*op for op in B.Ops])
end
####
struct ProductOfAbstractOperators{T} <: AbstractOperator where T <: AbstractOperator
    Ops::Vector{T}
end

struct SumOfAbstractOperators{T} <: AbstractOperator where T <: AbstractOperator
    Ops::Vector{T}
    c::Vector
end

function *(c::Number,L::Operator)
    SumOfAbstractOperators([L],[c])
end

function -(L::Operator)
    (-1)*L
end

function -(L::BlockAbstractOperator)
    BlockAbstractOperator(-L.Ops)
end

function +(Op1::AbstractOperator,Op2::AbstractOperator)
    SumOfAbstractOperators([Op1;Op2],[1;1])
end

function -(Op1::AbstractOperator,Op2::AbstractOperator)
    SumOfAbstractOperators([Op1;Op2],[1;-1])
end

function +(Op1::BlockAbstractOperator,Op2::BlockAbstractOperator)
    BlockAbstractOperator(Op1.Ops .+ Op2.Ops)
end

function -(Op1::BlockAbstractOperator,Op2::BlockAbstractOperator)
    BlockAbstractOperator(Op1.Ops .- Op2.Ops)
end

function +(S1::SumOfAbstractOperators{T1},S2::SumOfAbstractOperators{T2}) where {T1 <: AbstractOperator, T2 <:AbstractOperator}
    SumOfAbstractOperators(vcat(S1.Ops,S2.Ops),vcat(S1.c,S2.c))
end

function +(S1::AbstractOperator,S2::SumOfAbstractOperators{T2}) where {T2 <:AbstractOperator}
    (1*S1) + S2
end

function +(S2::SumOfAbstractOperators{T2},S1::AbstractOperator) where {T2 <:AbstractOperator}
    S2 + (1*S1)
end

struct AbstractZeroOperator <: AbstractOperator end
AbstractZeroOperator(n::Int64,m::Int64) = n == m && n == 1 ? AbstractZeroOperator() : BlockAbstractOperator(fill(AbstractZeroOperator(),n,m))
zero(::Type{AbstractOperator}) = AbstractZeroOperator()
zero(::AbstractOperator) = AbstractZeroOperator()

function *(Op::AbstractZeroOperator,Op2::AbstractOperator)
    Op
end

function +(Op::AbstractZeroOperator,Op2::AbstractOperator)
    Op2
end

function *(Op::AbstractOperator,Op2::AbstractZeroOperator)
    Op2
end

function +(Op::AbstractOperator,Op2::AbstractZeroOperator)
    Op
end

function +(Op::AbstractZeroOperator,Op2::AbstractZeroOperator)
    Op
end

function +(Op::SumOfAbstractOperators{T2}, Op2::AbstractZeroOperator) where T2<:AbstractOperator
    Op
end

struct Derivative <: AbstractOperator
    order::Integer
end

struct FloquetDerivative <: AbstractOperator
    order::Integer
    μ::Float64
end

struct Multiplication <: AbstractOperator
    f::Union{Function,BasisExpansion,Vector}
end

### SEMI ABSTRACT OPERATORS ###
struct Conversion <: AbstractOperator
    range::Basis
end

struct FastConversion <: AbstractOperator
    range::Basis
end

struct CoefConversion <: AbstractOperator
    range::Basis
end


struct BoundaryValue <: AbstractOperator
    o::Int64
    range::Basis
end
function BoundaryValue(o::Int64,ran::DirectSum)
    BlockDiagonalAbstractOperator(map(z -> BoundaryValue(o,z), ran.bases))
end

struct Residue <: AbstractOperator
    range::Basis
end
### Wrappers ###
struct Truncation <: AbstractOperator
    Op::AbstractOperator
    k::Int64
end
### ###

struct CauchyTransform <: AbstractOperator end

struct CauchyOperator <: AbstractOperator
    o::Int64
end

struct FourierTransform <: AbstractOperator
    o::Int64
end

Derivative() = Derivative(1)
FloquetDerivative(μ) = FloquetDerivative(1,μ)

function *(D1::Derivative,D2::Derivative)
    Derivative(D1.order + D2.order)
end

function ^(D1::Derivative,k::Integer)
    Derivative(k*D1.order)
end


function *(D1::FloquetDerivative,D2::FloquetDerivative)
    if D1.μ ≈ D2.μ
        return FloquetDerivative(D1.order + D2.order,D1.μ)
    else
        @error "Floquet parameters do not match"
    end
end

function ^(D1::FloquetDerivative,k::Integer)
    FloquetDerivative(k*D1.order,D1.μ)
end

function *(M::AbstractOperator,Op2::AbstractOperator)
    ProductOfAbstractOperators([M;Op2])
end

function *(P::ProductOfAbstractOperators,Op::AbstractOperator)
    ProductOfAbstractOperators(vcat(P.Ops,[Op]))
end

function *(Op::ProductOfAbstractOperators,sp::Basis)
    p = Op.Ops[end]*sp
    for i = length(Op.Ops)-1:-1:1
        p = Op.Ops[i]*p
    end
    p
end

function *(Op::SumOfAbstractOperators,sp::Basis)
    ops = [op*sp for op in Op.Ops]
    L = SumOfMatrixOperators([op.L for op in ops],Op.c)
    ConcreteOperator(ops[1].domain,ops[1].range,L)
end

function *(Op::AbstractOperator,f::BasisExpansion)
    Opc = Op*f.basis
    Opc*f
end

function *(Op::AbstractZeroOperator,b::Basis)
    ConcreteOperator(AnyBasis(),AnyBasis(),ZeroOperator())
end

function *(Op::AbstractZeroOperator,b::DirectSum)
    m = b.bases |> length
    B = BlockAbstractOperator(fill(AbstractZeroOperator(),m,m))
    B*b
end

function *(Op::AbstractOperator,b::DirectSum)
    B = [Op for bb in b.bases] |> diagm
    B*b
end

function *(Op::BlockAbstractOperator,sp::Basis)
    if size(Op.Ops)[2] > 1
        @error "Incorrect block size."
        return
    end
    COps = [op*sp for op in Op.Ops]
    sps = [op.range for op in COps][:,1]
    Ls = [op.L for op in COps]
    ConcreteOperator(sp,DirectSum(sps),BlockMatrixOperator(Ls))
end