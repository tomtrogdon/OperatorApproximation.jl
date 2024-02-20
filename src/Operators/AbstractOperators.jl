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

struct ProductOfAbstractOperators{T} <: AbstractOperator where T
    Ops::Vector{T}
end

struct SumOfAbstractOperators{T} <: AbstractOperator where T
    Ops::Vector{T}
    c::Vector
end

function *(c::Number,L::Operator)
    SumOfAbstractOperators([L],[c])
end

function -(L::Operator)
    (-1)*L
end

function +(Op1::AbstractOperator,Op2::AbstractOperator)
    SumOfAbstractOperators([Op1;Op2],[1;1])
end

function -(Op1::AbstractOperator,Op2::AbstractOperator)
    SumOfAbstractOperators([Op1;Op2],[1;-1])
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

struct Derivative <: AbstractOperator
    order::Integer
end

function Evaluation()
    Conversion(gridvals)
end

function Evaluation(b::T) where T <: GridValues
    Conversion(b)
end

struct Multiplication <: AbstractOperator
    f::Function
end

### SEMI ABSTRACT OPERATOR ###
struct Conversion <: AbstractOperator
    range::Basis
end
###

Derivative() = Derivative(1)

function *(D1::Derivative,D2::Derivative)
    Derivative(D1.order + D2.order)
end

function ^(D1::Derivative,k::Integer)
    Derivative(k*D1.order)
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
    SumOfConcreteOperators(ops[1].domain,ops[1].range,ops ,Op.c)
end

function *(Op::AbstractOperator,f::BasisExpansion)
    Opc = Op*f.basis
    Opc*f
end