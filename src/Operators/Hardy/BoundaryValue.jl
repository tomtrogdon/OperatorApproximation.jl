function *(B::BoundaryValue,b1::Hardy{T,S}) where {T <: Exterior, S <: Interval}
    if !(b1.GD.D == B.range.GD.D)
        @error "incompatible domain"
        return
    end

    if typeof(B.range.GD) <: DirectedLobattoMappedInterval
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
    ConcreteLazyOperator(b1,B.range,Op)
end