function isconvertible(b1::Hermite,b2::DiscreteBasis)
    iscompatible(b1.GD,b2.GD)
end

function conversion(b1::HermitePoly,b2::GridValues)
    basegrid =  n -> b2.GD.grid(n)
    # In principle, we need to do this:
    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    # but we are checking that the two grid domains are compatible
    # and currently this forces the composition of the maps to
    # be the identity
    a, b = Hermite_ab()
    Op = OPEvaluationOperator(basegrid,a,b)
    ConcreteLazyOperator(b1,b2,Op)
end

function conversion(b1::HermiteFun,b2::GridValues)
    basegrid =  n -> b2.GD.grid(n)
    # In principle, we need to do this:
    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    # but we are checking that the two grid domains are compatible
    # and currently this forces the composition of the maps to
    # be the identity
    a, b = Hermite_ab()
    Op = OPWeightedEvaluationOperator(basegrid,a,b,x -> (2pi)^(-0.25)exp(-x^2/4))
    ConcreteLazyOperator(b1,b2,Op)
end

function conversion(b1::HermitePoly,b2::FixedGridValues)
    # See conversion remark above.
    a, b = Hermite_ab()
    Op = FixedGridOPEvaluationOperator(b2.pts,a,b)
    ConcreteLazyOperator(b1,b2,Op)
end

function conversion(b1::HermiteFun,b2::FixedGridValues)
    # See conversion remark above.
    a, b = Hermite_ab()
    Op = FixedGridWeightedOPEvaluationOperator(b2.pts,a,x -> (2pi)^(-0.25)exp(-x^2/4))
    ConcreteLazyOperator(b1,b2,Op)
end