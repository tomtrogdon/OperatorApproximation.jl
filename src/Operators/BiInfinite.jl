function Matrix(Op::SemiLazyBandedOperator{T},n,m) where  T <: ZZ
    if n > size(Op.A)[1] || m > size(Op.A)[2]
        Op.A = Op.mat(max(n,m))
    end
    N,M = size(Op.A)
    Mm = N₋(M)
    mm = N₋(m)
    Nm = N₋(N)
    nm = N₋(n)
    Op.A(Nm - nm + 1: Nm - nm + 1 + n, Mm - mm + 1: Mm - mm + 1 + m)
end

function getlen(i,n,m) # suppose n <= m
    if 0 <= i <= m - n -1
        return n
    elseif i > m - n - 1
        return n - (i - (m - n - 1)) + 1
    else
        return n + i
    end
end

function getmat(F,n,m,nm,np)
    mM = OperatorApproximation.N₋(m)
    nM = OperatorApproximation.N₋(n)
    A = spzeros(n,m)
    if nM <= mM
        zero_diag = mM - nM
        for i = 0:m-1
            if -nm + zero_diag <= i <= zero_diag + np
                A += spdiagm(n,m, i => [F(-nM + j, -mM + j + i) for j = 0:getlen(i,n,m)-1])
            end
        end
        for i = -n+1:-1
            if -nm + zero_diag <= i <= zero_diag + np
                A += spdiagm(n,m, i => [F(-nM + j - i, -mM + j) for j = 0:getlen(i,n,m)-1])
            end
        end
        A
    else
        A = getmat((i,j) -> F(j,i),m,n) |> transpose
    end
    A
end

function Matrix(Op::BasicBandedOperator{T},n,m) where T <: ZZ
    mM = N₋(m)
    nM = N₋(n)
    A = spzeros(n,m)
    if nM <= mM
        zero_diag = mM - nM
        for i = 0:m-1
            if -Op.nm + zero_diag <= i <= zero_diag + Op.np
                A += spdiagm(n,m, i => [Op.A(-nM + j, -mM + j + i) for j = 0:getlen(i,n,m)-1])
            end
        end
        for i = -n+1:-1
            if -Op.nm + zero_diag <= i <= zero_diag + Op.np
                A += spdiagm(n,m, i => [Op.A(-nM + j - i, -mM + j) for j = 0:getlen(i,n,m)-1])
            end
        end
        A
    else
        A = getmat((i,j) -> Op.A(j,i),m,n) |> transpose
    end
    A
end

function Matrix(Op::MultipliedBandedOperator{T},n,m) where T <: ZZ  # TEST THIS
    cols = m
    ex = max(Op.V[end].nm,Op.V[end].np)
    rows = max(cols+2ex,1)
    A = Matrix(Op.V[end],rows,cols)
    for j = length(Op.V)-1:-1:2
        cols = rows
        ex = max(Op.V[j].nm,Op.V[j].np)
        rows = max(cols + 2ex,1)
        A = Matrix(Op.V[j],rows,cols)*A
    end
    cols = rows
    rows = n
    A = Matrix(Op.V[1],n,cols)*A
    A
end

## TODO: Multiplication routines could be possibly simplified
for op in (:BasicBandedOperator,:SemiLazyBandedOperator)
    @eval begin 
        function *(Op::$op{T},c::Vector) where T <: ZZ
            n = length(c)
            ex = max(Op.nm,Op.np)
            m = max(2ex + n,1)
            Matrix(Op,m,n)*c
        end

        function rowgrowth(Op::$op{T}) where T <: ZZ
            2*max(Op.nm,Op.np)
        end
    end
end

function *(Op::MultipliedBandedOperator{T},c::Vector) where T <: ZZ
    cols = length(c)
    ex = max(Op.V[end].nm,Op.V[end].np)
    rows = max(cols+2ex,1)
    v = Matrix(Op.V[end],rows,cols)*c
    for j = length(Op.V)-1:-1:1
        cols = rows
        ex = max(Op.V[j].nm,Op.V[j].np)
        rows = max(cols + 2ex,1)
        v = Matrix(Op.V[j],rows,cols)*v
    end
    v
end

function rowgrowth(Op::MultipliedBandedOperator{T}) where T <: ZZ
    cols = 0
    ex = max(Op.V[end].nm + Op.V[end].np)
    rows = cols+2ex
    for j = length(Op.V)-1:-1:1
        cols = rows
        ex = max(Op.V[j].nm + Op.V[j].np)
        rows = cols + 2ex
    end
    rows
end

function BIIdentityOperator()
    BasicBandedOperator(BI,0,0, (i,j) -> i == j ? 0.0 : 1)
end