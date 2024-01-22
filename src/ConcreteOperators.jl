abstract type LazyOperator <: Operator end

struct ConcreteOperator{D<:Basis,R<:Basis} <: Operator
    domain::D
    range::R
    L::LazyOperator
end

struct SumOfConcreteOperators{T} <: Operator where T
    Ops::Vector{T}
    c::Vector
end

struct BandedOperator <: LazyOperator
    nm::Integer
    np::Integer
    A::Function
end

struct GridMultiplication <: LazyOperator
    GV::GridValues
    f::Function
end

abstract type BasisEvaluationOperator <: LazyOperator end

struct OPEvaluationOperator <: BasisEvaluationOperator  ## Add CollocatedOperator?
    grid::Function
    a::Function
    b::Function
end

function +(Op1::ConcreteOperator,Op2::ConcreteOperator) # need to check that the range and domain are compatible.
    SumOfConcreteOperators([Op1;Op2],[1;1])
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

struct VariableCollocatedBandedOperator <: LazyOperator
    M::GridMultiplication
    Op::LazyOperator
end

function *(M::GridMultiplication,Op::CollocatedBandedOperator)
    VariableCollocatedBandedOperator(M,Op)
end

function *(M::GridMultiplication,Op::OPEvaluationOperator)
    VariableCollocatedBandedOperator(M,Op)
end

function *(VC::VariableCollocatedBandedOperator,L::LazyOperator)
    VariableCollocatedBandedOperator(VC.M,VC.Op*L)
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

function *(M::CollocatedMultiplication, Op::ConcreteOperator{D,R}) where {D, R <: GridValues}
    MM = GridMultiplication(Op.range,M.f)
    ConcreteOperator(Op.domain,Op.range,MM*Op.L)
end

function *(M::Multiplication, Op::ConcreteOperator{D,R}) where {D, R <: GridValues}
    MM = GridMultiplication(Op.range,M.f)
    ConcreteOperator(Op.domain,Op.range,MM*Op.L)
end

function Matrix(Op::VariableCollocatedBandedOperator,n,m)  # need to map the values.
    Diagonal(Op.M.GV.GD.grid(n))*Matrix(Op.Op,n,m)
end

function Matrix(Op::SumOfConcreteOperators,n,m)
    sum( [Op.c[i]*Matrix(Op.Ops[i],n,m) for i in length(Op.c)])
end