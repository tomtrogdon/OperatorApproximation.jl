function toeplitz(c::Vector,N::Integer)
    mm = length(c) # should be simpler
    m₋ = N₋(mm)
    m₊ = N₊(mm)
    dm = N - m₋ - 1
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

function toeplitz_function(c::Vector)
    function F(k,j)
        n = length(c)
        np = N₊(length(c))
        nm = N₋(length(c))
        # F(0,0) = c[nm+1]
        # F(j,j) = c[nm+1]
        # F(j - i, j) = c[nm + 1 - i], k = j - i, i = j - k
        # F(k,j) = c[nm + 1 + j -k]
        if np + 1 + j -k< 1 || np + 1 + j -k > length(c)
            return 0.0
        else
            return c[nm - j + k + 1]
        end
    end
end


function *(M::FastMultiplication,sp::Fourier)
    if typeof(M.f) <: Function
        a = sp.GD.D.a
        b = sp.GD.D.b
        GD = PeriodicMappedInterval(a,b)
        ff = BasisExpansion(M.f,Fourier(GD)) |> chop
    else
        ff = M.f
    end
    growth = M.growth < 0 ? length(ff.c) : M.growth
    f_grid = n -> midft(pad(ℤ,ff.c,n))
    Op = FastGridMultiplication(f_grid, growth)
    ConcreteOperator(sp,sp,Op)
end

function *(M::Multiplication,sp::Fourier)
    if typeof(M.f) <: Function
        a = sp.GD.D.a
        b = sp.GD.D.b
        GD = PeriodicMappedInterval(a,b)
        ff = BasisExpansion(M.f,Fourier(GD)) |> chop
    else 
        ff = M.f
    end
    
    if typeof(ff.basis) <: Fourier && isconvertible(ff.basis,sp)
        np = N₋(length(ff.c)); nm = length(ff.c) - np + 1
        Op = BasicBandedOperator{ℤ,ℤ}(np,np,toeplitz_function(ff.c))
    else 
        1 + 1 #TODO: just evaluate and expand, need transform
    end
    ConcreteOperator(sp,sp,Op)
end