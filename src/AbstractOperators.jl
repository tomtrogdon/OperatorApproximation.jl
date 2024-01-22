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

function +(Op1::AbstractOperator,Op2::AbstractOperator)
    SumOfAbstractOperators([Op1;Op2],[1;1])
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

struct Evaluation <: AbstractOperator
    GD::GridDomain
end

Evaluation() = Evaluation(NullGrid())

ℰ = x -> Evaluation(x)
𝒟 = x -> Derivative(x)

struct Multiplication <: AbstractOperator
    f::Function
end

ℳ = f -> Multiplication(f)

struct CollocatedOperator <: AbstractOperator
   Op::AbstractOperator
   GD::GridDomain
end

struct CollocatedMultiplication <: AbstractOperator
    f::Function
end

Derivative() = Derivative(1)

function *(D1::Derivative,D2::Derivative)
    Derivative(D1.order + D2.order)
end

function *(E::Evaluation,Op::ProductOfAbstractOperators)
    if typeof(Op.Ops[1]) <: Multiplication
        PoO = ProductOfAbstractOperators(Op.Ops[2:end])
        CollocatedMultiplication(Op.Ops[1].f)*(E*PoO)
    else
        CollocatedOperator(Op,E.GD)
    end
end

function *(E::Evaluation,Op::AbstractOperator)
    CollocatedOperator(Op,E.GD)
end

function *(M::Multiplication,Op2::AbstractOperator)
    ProductOfAbstractOperators([M;Op2])
end

function *(M::Multiplication,Op::CollocatedOperator)
    ProductOfAbstractOperators([CollocatedMultiplication(M.f);Op])
end

function *(Op1::CollocatedOperator,Op2::AbstractOperator)
    CollocatedOperator(Op1.Op*Op2,Op1.GD)
end

function *(E::Evaluation,M::Multiplication)
    ProductOfAbstractOperators([M;E])  # Note that this is an approximation.
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
    SumOfConcreteOperators([op*sp for op in Op.Ops],Op.c)
end
