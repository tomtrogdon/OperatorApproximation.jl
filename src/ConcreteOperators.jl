abstract type LazyOperator <: Operator end

struct ConcreteOperator{D<:Basis,R<:Basis} <: Operator
    domain::D
    range::R
    L::LazyOperator
end

struct SumOfConcreteOperators{D<:Basis,R<:Basis,T} <: Operator where T
    domain::D
    range::R
    Ops::Vector{T}
    c::Vector
end

abstract type BandedOperator <: LazyOperator end

struct SingleBandedOperator <: LazyOperator
    nm::Integer
    np::Integer
    A::Function
end

abstract type GridOperator <: LazyOperator end

# struct DiagonalLazyOperator <: GridOperator
#     diag::Function
# end

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

function +(Op1::ConcreteOperator,Op2::ConcreteOperator) # need to check that the range and domain are compatible.
    #Check domain and range here    
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

# function ToDiag(M::GridMultiplication)
#     DiagonalLazyOperator(n -> M.f.(M.GV.GD.grid(n)))
# end

function *(M::GridMultiplication,Op::CollocatedBandedOperator)
    GridTimesCollocated([M],Op)
end

function *(M::GridMultiplication,Op::LazyCollocatedOperator)
    GridTimesCollocated([M],Op)
end

# function *(M::GridMultiplication,Op::OPEvaluationOperator)
#     ToDiag(M)*Op
# end

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

function Matrix(Op::SumOfConcreteOperators,n,m)
    # TODO: check domain & range
    sum( [Op.c[i]*Matrix(Op.Ops[i],n,m) for i in length(Op.c)])
end

function changegrid(C::ConcreteOperator{GridValues,Ran},GV::GridValues) where Ran <: Basis
    ConcreteOperator(C.domain,GV,C.L)
end

function changegrid(C::SumOfConcreteOperators{GridValues,Ran},GV::GridValues) where Ran <: Basis
    ConcreteOperator(C.domain,GV,C.L)
end