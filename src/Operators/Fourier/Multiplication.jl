function toeplitz(c::Vector,N::Integer)
    mm = length(c) # should be simpler
    m₋ = N₋(mm)
    m₊ = N₊(mm)
    dm = N - m₋ - 1
    display(dm)
    range = m₋:-1:-convert(Int64,floor((mm-1)/2))
    range = range .+ (dm < 0 ? dm : 0)
    mstart = max(0,-dm)
    mend = min(m₋ + dm,m₊)
    diags1 = [fill(c[1+j],1+j+dm) for j in mstart:m₋]
    diags2 = [fill(c[m₋+1+j],m₋+1-j+dm) for j in 1:mend]
    diags = vcat(diags1,diags2)
    inds = map( (x,y) -> x => y, range,diags)
    A = spdiagm(inds[1:1]...)
    for i=2:length(inds)
        A += spdiagm(inds[i:i]...)
    end
    A
end




function *(M::Multiplication,sp::Fourier)
    if typeof(M.f) <: Function
        a = sp.GD.D.a
        b = sp.GD.D.b
        GD = PeriodicMappedInterval(a,b,0.0)
        ff = BasisExpansion(M.f,Fourier(GD),200) |> Chop
    else 
        ff = M.f
    end
    
    if typeof(ff.basis) <: Fourier && isconvertible(ff.basis,sp)
        

        Op = LazyBandedOperator(m,m,f,f(5))
    else 
        1 + 1 #TODO: just evaluate and expand, need transform
    end
    ConcreteLazyOperator(sp,sp,Op)
end