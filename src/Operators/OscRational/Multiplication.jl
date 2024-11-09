function *(M::Multiplication,sp::OscRational)
    if typeof(M.f) <: Function 
        GD = RationalRealAxis()
        α = 0.0
        ff = BasisExpansion(M.f,OscRational(GD,α)) |> chop
    else 
        ff = M.f
    end
    if typeof(ff.basis) <: OscRational && isconvertible(ff.basis,sp)
        np = N₋(length(ff.c)); nm = length(ff.c) - np + 1 #why even define nm??
        Op = BasicBandedOperator{ℤ,ℤ}(np,np,toeplitz_function(ff.c)) #creates Toeplitz operator
    else 
        1 + 1 #TODO: just evaluate and expand, need transform #I am assumping this will just use toeplitz()
    end

    sp2 = OscRational(RationalRealAxis(),α + sp.α)
    ConcreteOperator(sp,sp2,Op) #operator in practice for multiplication
end