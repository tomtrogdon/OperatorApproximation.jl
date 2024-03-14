function *(M::Multiplication,sp::Laurent)
    if typeof(M.f) <: Function
        cen = sp.GD.D.cen
        rad = sp.GD.D.rad
        GD = PeriodicMappedCircle(cen,rad)
        ff = BasisExpansion(M.f,Laurent(GD)) |> chop
    else 
        ff = M.f
    end
    
    if typeof(ff.basis) <: Laurent && isconvertible(ff.basis,sp)
        np = N₋(length(ff.c)); nm = length(ff.c) - np + 1
        Op = BasicBandedOperator(BI,np,np,toeplitz_function(ff.c))
    else 
        1 + 1 #TODO: just evaluate and expand, need transform
    end
    ConcreteLazyOperator(sp,sp,Op)
end