function isconvertible(b1::Ultraspherical,b2::DiscreteBasis)
    iscompatible(b1.GD,b2.GD)
end

function isconvertible(b1::Ultraspherical,b2::Ultraspherical)
    iscompatible(b1.GD,b2.GD)
end

function conversion(b1::Ultraspherical,b2::GridValues)
    basegrid =  n -> b2.GD.grid(n)
    # In principle, we need to do this:
    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    # but we are checking that the two grid domains are compatible
    # and currently this forces the composition of the maps to
    # be the identity
    a, b = Jacobi_ab(b1.λ - 1/2, b1.λ - 1/2)
    Op = OPEvaluationOperator(basegrid,a,b)
    ConcreteLazyOperator(b1,b2,Op)
end

function conversion(b1::Ultraspherical,b2::FixedGridValues)
    # See conversion remark above.
    a, b = Jacobi_ab(b1.λ - 1/2, b1.λ - 1/2)
    Op = FixedGridOPEvaluationOperator(b2.pts,a,b)
    ConcreteLazyOperator(b1,b2,Op)
end