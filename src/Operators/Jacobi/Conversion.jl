function isconvertible(b1::Jacobi,b2::DiscreteBasis)
    iscompatible(b1.GD,b2.GD)
end

function conversion(b1::Jacobi,b2::GridValues)
    basegrid =  n -> b2.GD.grid(n)
    # In principle, we need to do this:
    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    # but we are checking that the two grid domains are compatible
    # and currently this forces the composition of the maps to
    # be the identity
    a, b = Jacobi_ab(b1.α, b1.β)
    Op = OPEvaluationOperator(basegrid,a,b)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Jacobi,b2::FixedGridValues)
    # See conversion remark above.
    a, b = Jacobi_ab(b1.α, b1.β)
    Op = OPEvaluationOperator(b2.pts,a,b)
    ConcreteOperator(b1,b2,Op)
end