abstract type Operator end

abstract type AbstractOperator <: Operator end

struct SumOfOperators{T} <: Operator where T
    Ops::Vector{T}
    c::Vector
end

function *(c::Number,L::Operator)
    SumOfOperators([L],[c])
end

function +(S1::SumOfOperators{T1},S2::SumOfOperators{T2}) where {T1 <: AbstractOperator, T2 <:AbstractOperator}
    SumOfOperators(vcat(S1.Ops,S2.Ops),vcat(S1.c,S2.c))
end

function +(S1::AbstractOperator,S2::SumOfOperators{T2}) where {T2 <:AbstractOperator}
    (1*S1) + S2
end

function +(S2::SumOfOperators{T2},S1::AbstractOperator) where {T2 <:AbstractOperator}
    S2 + (1*S1)
end

struct Derivative <: AbstractOperator
    order::Integer
end

struct Evaluation <: AbstractOperator end

struct Multiplication <: AbstractOperator
    f::Function
end

struct CollocatedOperator <: AbstractOperator
   Op::AbstractOperator
end

Derivative() = Derivative(1)

function *(D1::Derivative,D2::Derivative)
    Derivative(D1.order + D2.order)
end

function *(E::Evaluation,Op::AbstractOperator)
    CollocatedOperator(Op)
end

function *(Op1::CollocatedOperator,Op2::AbstractOperator)
    CollocatedOperator(Op1.Op*Op2)
end

abstract type LazyOperator <: Operator end

struct ConcreteOperator{D<:Basis,R<:Basis} <: Operator
    domain::D
    range::R
    L::LazyOperator
end

struct BandedOperator <: LazyOperator
    nm::Integer
    np::Integer
    A::Function
end

abstract type BasisEvaluationOperator <: LazyOperator end

struct OPEvaluationOperator <: BasisEvaluationOperator  ## Add CollocatedOperator?
    grid::Function
    a::Function
    b::Function
end

function Matrix(Op::OPEvaluationOperator,n,m)
    poly(Op.a,Op.b,m,Op.grid(n)) 
 end

 function Matrix(Op::BandedOperator,n,m)
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

struct MultipliedBandedOperator <: LazyOperator
    V::Vector{BandedOperator}
end

function *(A::BandedOperator,B::BandedOperator)
    MultipliedBandedOperator([A;B])
end

function *(A::BandedOperator,B::MultipliedBandedOperator)
    MultipliedBandedOperator(vcat([A],B.V))
end
    
function *(B::MultipliedBandedOperator,A::BandedOperator)
    MultipliedBandedOperator(vcat(B.V,[A]))
end
        
function *(B::MultipliedBandedOperator,A::MultipliedBandedOperator)
    MultipliedBandedOperator(vcat(B.V,A.V))
end

struct CollocatedBandedOperator <: LazyOperator
    V::Vector{BandedOperator}
    E::BasisEvaluationOperator
end

function *(E::BasisEvaluationOperator,A::BandedOperator)
    CollocatedBandedOperator([A],E)
end

function *(E::BasisEvaluationOperator,A::MultipliedBandedOperator)
    CollocatedBandedOperator(A.V,E)
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

function Matrix(Op::ConcreteOperator,n,m)
    Matrix(Op.L,n,m)
end

function Matrix(Op::CollocatedBandedOperator,n,m)
    cols = m
    rows = max(cols+Op.V[end].nm,1)
    A = Matrix(Op.V[end],rows,cols)
    for j = length(Op.V)-1:-1:1
        cols = rows
        rows = max(cols + Op.V[j].nm,1)
        A = Matrix(Op.V[j],rows,cols)*A
    end
    cols = rows
    rows = n
    A = Matrix(Op.E,rows,cols)*A
    A
end