function *(M::Multiplication,sp::OscRational)
    if typeof(M.f) <: Function 
        GD = sp.GD
        α = 0.0
        ff = BasisExpansion(M.f,OscRational(GD,α)) |> chop
    else 
        ff = M.f
        α = ff.basis.α
    end

    if typeof(ff.basis) <: OscRational && iscompatible(ff.basis.GD,sp.GD)
        np = N₋(length(ff.c)); nm = length(ff.c) - np + 1 #why even define nm??
        Op = BasicBandedOperator{ℤ,ℤ}(np,np,toeplitz_function(ff.c)) #creates Toeplitz operator
    else 
        1 + 1 #TODO: just evaluate and expand, need transform #I am assumping this will just use toeplitz()
    end

    if α ≈ 0.0
        return ConcreteOperator(sp,sp,Op)
    else
        sp2 = OscRational(sp.GD,α + sp.α)
        return ConcreteOperator(sp,sp2,Op) #operator in practice for multiplication
    end
end