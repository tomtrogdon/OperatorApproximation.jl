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

function *(Op::SumOfAbstractOperators,C::ConcreteLazyOperator)
    ops = [op*C for op in Op.Ops]
    SumOfConcreteOperators(ops[1].domain,ops[1].range,ops,Op.c)
end

abstract type BandedOperator <: LazyOperator end

struct SingleBandedOperator <: LazyOperator
    nm::Integer
    np::Integer
    A::Function
end

abstract type GridOperator <: LazyOperator end

struct GridMultiplication <: GridOperator
    GV::GridValues
    f::Function
end

abstract type LazyCollocatedOperator <: LazyOperator end
abstract type BasisEvaluationOperator <: LazyCollocatedOperator end

struct GridTimesCollocated <: LazyOperator
    GOps::Vector{GridOperator}
    C::LazyCollocatedOperator
end

struct OPEvaluationOperator <: BasisEvaluationOperator  ## Add CollocatedOperator?
    grid::Function
    a::Function
    b::Function
end

struct ConcreteProjector{T} <: Operator where T <: Basis
    domain::T
    N::Integer
end

function *(P::Projector,basis::Basis)
    ConcreteProjector(basis,P.N)
end

function (P::Projector)(Op::ConcreteOperator)
    ConcreteProjector(Op.range,P.N)
end

function *(P::ConcreteProjector{T},f::Function) where T <: GridValues
    grid = P.domain.GD.grid(P.N)
    grid = P.domain.GD.D.map.(grid)
    f.(grid)
end

abstract type ConcreteFunctional <: ConcreteOperator end

struct ConcreteBoundaryFunctional{T<:Basis} <: ConcreteFunctional
    domain::T
    A::Matrix
    B::Matrix
end

function *(BF::BoundaryFunctional,sp::Basis)
    ConcreteBoundaryFunctional(sp,BF.A,BF.B)
end

function +(Op1::ConcreteLazyOperator,Op2::ConcreteLazyOperator) # need to check that the range and domain are compatible.
    #TODO: Check domain and range here    
    SumOfConcreteOperators(Op1.domain,Op1.range,[Op1;Op2],[1;1])
end

function Matrix(Op::OPEvaluationOperator,n,m)
    poly(Op.a,Op.b,m,Op.grid(n)) 
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

struct CollocatedBandedOperator <: LazyCollocatedOperator
    V::Vector{SingleBandedOperator}
    E::BasisEvaluationOperator
end

function *(E::BasisEvaluationOperator,A::SingleBandedOperator)
    CollocatedBandedOperator([A],E)
end

function *(E::BasisEvaluationOperator,A::MultipliedBandedOperator)
    CollocatedBandedOperator(A.V,E)
end

function *(M::GridMultiplication,Op::CollocatedBandedOperator)
    GridTimesCollocated([M],Op)
end

function *(M::GridMultiplication,Op::LazyCollocatedOperator)
    GridTimesCollocated([M],Op)
end

function Matrix(M::GridMultiplication,n,m)
    nm = max(n,m)
    X = M.GV.GD.D.map.(M.GV.GD.grid(nm))
    A = Diagonal(M.f.(X))
    A[1:n,1:m] |> sparse
end

function Matrix(Op::GridTimesCollocated,n,m)
    A = Matrix(Op.C,n,m)
    for i = length(Op.GOps):-1:1
        A = Matrix(Op.GOps[i],n,n)*A
    end
    A
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

function Matrix(Op::ConcreteLazyOperator,n,m)
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

function *(M::CollocatedMultiplication, Op::ConcreteLazyOperator{D,R}) where {D, R <: GridValues}
    MM = GridMultiplication(Op.range,M.f)
    ConcreteLazyOperator(Op.domain,Op.range,MM*Op.L)
end

function *(M::Multiplication, Op::ConcreteLazyOperator{D,R}) where {D, R <: GridValues}
    MM = GridMultiplication(Op.range,M.f)
    ConcreteLazyOperator(Op.domain,Op.range,MM*Op.L)
end

function Matrix(Op::SumOfConcreteOperators,n,m)
    # TODO: check domain & range
    A = Op.c[1]*Matrix(Op.Ops[1],n,m)
    for i = 2:length(Op.c)
        A += Op.c[i]*Matrix(Op.Ops[i],n,m)
    end
    A
end

function changegrid(C::ConcreteLazyOperator{GridValues,Ran},GV::GridValues) where Ran <: Basis
    ConcreteLazyOperator(C.domain,GV,C.L)
end

function changegrid(C::SumOfConcreteOperators{GridValues,Ran},GV::GridValues) where Ran <: Basis
    ConcreteLazyOperator(C.domain,GV,C.L)
end

function canonicalBC(kl,kr)
    k = kl + kr
    A = spzeros(k,k)
    B = spzeros(k,k)
    for i = 1:kl
        A[i,i] = 1.0
    end
    for i = 1:kr
        B[i+kr,i] = 1.0
    end
    (A,B)
end