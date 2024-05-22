#I think this just says we can go from a Rational to a DiscreteBasis
function isconvertible(b1::Rational,b2::DiscreteBasis)  # probably can be simplified
    iscompatible(b1.GD,b2.GD)
end

#I think this just says we can go from a Rational basis to another Rational basis
function isconvertible(b1::Rational,b2::Rational)
    iscompatible(b1.GD,b2.GD)
end

#This converts the Rational basis to values on the grid; note that GridValues is a DiscreteBasis
function conversion(b1::Rational,b2::GridValues)
    basegrid =  n -> b2.GD.grid(n)
    # In principle, we need to do this:
    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    # but we are checking that the two grid domains are compatible
    # and currently this forces the composition of the maps to
    # be the identity
    
    #!!!!!!!!!!!!!!Need to change to some kind of RationalEvaluationOperator???
    Op = FourierEvaluationOperator(basegrid) #QUESTION: What does "FourierEvalutationOperator" do and where defined?
    ConcreteOperator(b1,b2,Op)
end

#This does the same as above but for fixed values on the grid
function conversion(b1::Rational,b2::FixedGridValues)
    # See conversion remark above.
    Op = FixedGridFourierEvaluationOperator(b2.pts) #QUESTION: Where does "FixedGridFourierEvaluationOperator" live?
    ConcreteOperator(b1,b2,Op)
end

#This converts the Fourier basis to itself
function conversion(b1::Rational,b2::Rational)
    # TODO:  identity operator
    ConcreteOperator(b1,b2,IdentityOperator())
end
