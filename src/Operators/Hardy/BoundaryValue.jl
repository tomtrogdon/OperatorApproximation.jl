function *(B::BoundaryValue,b1::Hardy{Exterior{T},S}) where {T <: Union{JacobiMappedInterval,JacobiInterval}, S <: Interval}
    if !(b1.GD.D == B.range.GD.D)
        return Conversion(B.range)*b1
    end

    if typeof(B.range.GD) <: DirectedGridInterval
        basegrid =  n -> B.range.GD.dgrid(n)
    else
        basegrid =  n -> B.range.GD.grid(n)
    end

    α = b1.GD.GD.α
    β = b1.GD.GD.β
    a, b = Jacobi_ab(α,β)
    if B.o == 1
        Op = OPCauchyEvaluationOperator(basegrid, a, b, JacobiSeedPos(α,β))
    elseif B.o == -1
        Op = OPCauchyEvaluationOperator(basegrid, a, b, JacobiSeedNeg(α,β))
    end
    ConcreteOperator(b1,B.range,Op)
end

function *(B::BoundaryValue,b1::Hardy{Exterior{T},S}) where {T <: Union{MarchenkoPasturMappedInterval,MarchenkoPasturInterval}, S <: Interval}
    if !(b1.GD.D == B.range.GD.D)
        return Conversion(B.range)*b1
    end

    if typeof(B.range.GD) <: DirectedGridInterval
        basegrid =  n -> B.range.GD.dgrid(n)
    else
        basegrid =  n -> B.range.GD.grid(n)
    end

    α = b1.GD.GD.d
    a, b = MP_ab(d)
    if B.o == 1
        Op = OPCauchyEvaluationOperator(basegrid, a, b, MPSeedPos(d))
    elseif B.o == -1
        Op = OPCauchyEvaluationOperator(basegrid, a, b, MPSeedNeg(d))
    end
    ConcreteOperator(b1,B.range,Op)
end

function *(B::BoundaryValue,b1::Hardy{T,S}) where {T <: Exterior, S <: DiscreteDomain}
    if typeof(B.range) <: FixedGridValues # don't map
        Op = PoleResCauchyEvaluationOperator(B.range.GD.grid,b1.GD.grid)
        return ConcreteOperator(b1,B.range,Op)
    else
        basegrid =  n -> B.range.GD.grid(n)
        gridfun = n -> B.range.GD.D.map(basegrid(n))  # map
        Op = PoleResCauchyEvaluationOperator(gridfun,b1.GD.grid)
        return ConcreteOperator(b1,B.range,Op)
    end
end

function *(B::BoundaryValue,b1::Hardy{Interior{T},S}) where {T <: Union{PeriodicCircle,PeriodicMappedCircle}, S <: Circle}
    basegrid = B.range.GD.grid
    if B.o == 1
        Op = PosLaurentEvaluationOperator(basegrid)
    elseif B.o == -1
        Op = ZeroOperator()
    end
    ConcreteOperator(b1,B.range,Op)
end

function *(B::BoundaryValue,b1::Hardy{Exterior{T},S}) where {T <: Union{PeriodicCircle,PeriodicMappedCircle}, S <: Circle}
    basegrid = B.range.GD.grid
    if B.o == -1
        Op = NegLaurentEvaluationOperator(basegrid)
    elseif B.o == 1
        Op = ZeroOperator()
    end
    ConcreteOperator(b1,B.range,Op)
end