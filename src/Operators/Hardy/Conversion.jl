function isconvertible(b1::Hardy,b2::Basis)
    true # maybe some restrinctions on boundary intersections?
end

function conversion(b1::Hardy{Exterior{T},S},b2::GridValues) where {T <: Union{JacobiMappedInterval,JacobiInterval}, S <: Interval}
    if typeof(b2.GD) <: DirectedGridInterval
        basegrid =  n -> b2.GD.dgrid(n)
    else
        basegrid =  n -> b2.GD.grid(n)
    end
    # We need to do this:
    gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    α = b1.GD.GD.α
    β = b1.GD.GD.β
    a, b = Jacobi_ab(α,β)
    Op = OPCauchyEvaluationOperator(gridfun, a, b, JacobiSeed(α,β))
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Hardy{Exterior{T},S},b2::GridValues) where {T <: Union{MarchenkoPasturMappedInterval,MarchenkoPasturInterval}, S <: Interval}
    if typeof(b2.GD) <: DirectedGridInterval
        basegrid =  n -> b2.GD.dgrid(n)
    else
        basegrid =  n -> b2.GD.grid(n)
    end
    # We need to do this:
    gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    d = b1.GD.GD.d
    a, b = MP_ab(d)
    Op = OPCauchyEvaluationOperator(gridfun, a, b, MPSeed(d))
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Hardy{Exterior{T},S},b2::FixedGridValues) where {T <: Union{JacobiMappedInterval,JacobiInterval}, S <: Interval}
    basegrid = b2.GD.grid
    # We need to do this:
    mapgrid = b1.GD.D.imap(b2.GD.D.map(basegrid))
    α = b1.GD.GD.α
    β = b1.GD.GD.β
    a, b = Jacobi_ab(α,β)
    Op = OPCauchyEvaluationOperator(mapgrid, a, b, JacobiSeed(α,β))
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Hardy{Exterior{T},S},b2::FixedGridValues) where {T <: Union{MarchenkoPasturMappedInterval,MarchenkoPasturInterval}, S <: Interval}
    basegrid = b2.GD.grid
    # We need to do this:
    mapgrid = b1.GD.D.imap(b2.GD.D.map(basegrid))
    d = b1.GD.GD.d
    a, b = MP_ab(d)
    Op = OPCauchyEvaluationOperator(mapgrid, a, b, MPSeed(d))
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Hardy{T,S},b2::GridValues) where {T <: Exterior, S <: DiscreteDomain}
    basegrid =  n -> b2.GD.grid(n)
    gridfun = n -> b2.GD.D.map(basegrid(n))
    Op = PoleResCauchyEvaluationOperator(gridfun,b1.GD.grid)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Hardy{T,S},b2::FixedGridValues) where {T <: Exterior, S <: DiscreteDomain}
    Op = PoleResCauchyEvaluationOperator(b2.GD.grid,b1.GD.grid)
    ConcreteOperator(b1,b2,Op)
end
