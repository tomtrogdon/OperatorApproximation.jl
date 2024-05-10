function isconvertible(b1::Hardy,b2::Basis)
    true # maybe some restrinctions on boundary intersections?
end

function conversion(b1::Hardy{T,S},b2::GridValues) where {T <: Exterior, S <: Interval}
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
    ConcreteLazyOperator(b1,b2,Op)
end

function conversion(b1::Hardy{T,S},b2::GridValues) where {T <: Exterior, S <: DiscreteDomain}
    Op = PoleResCauchyEvaluationOperator(b2.GD.grid,b1.GD.grid)
    ConcreteLazyOperator(b1,b2,Op)
end