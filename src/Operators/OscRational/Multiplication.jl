function *(M::Multiplication,sp::OscRational)
    if typeof(M.f) <: Function 
        GD = RationalRealAxis()
        ff = BasisExpansion(M.f,OscRational(GD)) |> chop
    else 
        ff = M.f
    end
    if typeof(ff.basis) <: OscRational && isconvertible(ff.basis,sp)
        np = N₋(length(ff.c)); nm = length(ff.c) - np + 1 #why even define nm??
        Op = BasicBandedOperator{ℤ,ℤ}(np,np,toeplitz_function(ff.c)) #creates Toeplitz operator
    else 
        1 + 1 #TODO: just evaluate and expand, need transform #I am assumping this will just use toeplitz()
    end
    ConcreteOperator(sp,sp,Op) #operator in practice for multiplication
end