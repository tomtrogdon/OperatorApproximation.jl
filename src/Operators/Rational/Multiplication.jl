#I believe this builts a Toeplitz multiplication matrix with the coeffs out to N coeffs on either diag?
function toeplitz(c::Vector,N::Integer) #QUESTION: is this N even needed???
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
    T = spdiagm(inds[1:1]...)
    for i=2:length(inds)
        T += spdiagm(inds[i:i]...)
    end
    
    M = T .- (c*(transpose(ones(mm)))) .- ((sum(c))*I)
end

#It seems like it is returning the elements in a toeplitz matrix??
function toeplitz_function(c::Vector)
    function F(k,j)
        nm = N₋(length(c))
        d = sum(c)
        # F(0,0) = c[nm+1]
        # F(j,j) = c[nm+1]
        # F(j - i, j) = c[nm + 1 - i], k = j - i, i = j - k
        # F(k,j) = c[nm + 1 + j -k] (assuming k is row and j is column)
        if nm + 1 + j -k< 1 || nm + 1 + j -k > length(c)
            return 0.0 - c[k] - d
        else
            return c[nm + 1 + j - k] - c[k] - d
        end
    end
end

#Actually defines how to apply multiplication to the Rational space
function *(M::Multiplication,sp::Rational)
    if typeof(M.f) <: Function #QUESTION: What is M.f???
        a = sp.GD.D.a
        b = sp.GD.D.b
        GD = PeriodicMappedInterval(a,b) #QUESTION: Is this the one we want??
        ff = BasisExpansion(M.f,Rational(GD)) |> chop
    else 
        ff = M.f
    end
    
    if typeof(ff.basis) <: Rational && isconvertible(ff.basis,sp)
        np = N₋(length(ff.c)); nm = length(ff.c) - np + 1 #why even define nm??
        Op = BasicBandedOperator{ℤ,ℤ}(np,np,toeplitz_function(ff.c)) #creates Toeplitz operator
    else 
        1 + 1 #TODO: just evaluate and expand, need transform #I am assumping this will just use toeplitz()
    end
    ConcreteOperator(sp,sp,Op) #operator in practice for multiplication
end