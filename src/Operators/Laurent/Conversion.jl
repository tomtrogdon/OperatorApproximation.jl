function isconvertible(b1::Laurent,b2::DiscreteBasis)  # probably can be simplified
    iscompatible(b1.GD,b2.GD)
end

function isconvertible(b1::Laurent,b2::Laurent)
    iscompatible(b1.GD,b2.GD)
end

function conversion(b1::Laurent,b2::GridValues)
    basegrid =  n -> b2.GD.grid(n)
    # In principle, we need to do this:
    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    # but we are checking that the two grid domains are compatible
    # and currently this forces the composition of the maps to
    # be the identity
    Op = LaurentEvaluationOperator(basegrid)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Laurent,b2::FixedGridValues)
    # See conversion remark above.
    Op = FixedGridLaurentEvaluationOperator(b2.pts)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Laurent,b2::Laurent)
    # TODO:  identity operator
    # ConcreteOperator(b1,b2,IdentityOperator())
    ConcreteOperator(b1,b2,BasicBandedOperator{ℤ,ℤ}(0,0, (i,j) -> i == j ? complex(1.0) : 0.0im ))
end