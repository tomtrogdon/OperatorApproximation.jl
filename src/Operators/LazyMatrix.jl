function Matrix(Op::ZeroOperator,n,m)
    spzeros(n,m)
end

for op in (:BasicBandedOperator,:SemiLazyBandedOperator)
    @eval begin 
        rowgrowth(Op::$op{T, S}) where  {T <: ℕ₊, S <: ℕ₊} = Op.nm
        rowgrowth(Op::$op{T, S}) where  {T <: ℕ₋, S <: ℕ₋} = Op.np
        rowgrowth(Op::$op{T, S}) where  {T <: ℤ, S <: ℤ} = 2*max(Op.nm,Op.np)
        rowgrowth(Op::$op{T, S}) where  {T <: ℤ, S <: ℕ₊} = Op.nm
        rowgrowth(Op::$op{T, S}) where  {T <: ℤ, S <: ℕ₋} = Op.np
        rowgrowth(Op::$op{T, S}) where  {T <: ℕ₋, S <: ℤ} = 2*Op.np
        rowgrowth(Op::$op{T, S}) where  {T <: ℕ₊, S <: ℤ} = 2*Op.nm

        function *(Op::$op{T,S},c::Vector) where {T, S}
            n = length(c)
            m = max(Op.nm + n,1)
            Matrix(Op,m,n)*c
        end
    end
end

function rowgrowth(Op::ProductOfBandedOperators)
    [rowgrowth(op) for op in Op.V] |> sum
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
    mM = N₋(m)
    nM = N₋(n)
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
        A = getmat((i,j) -> F(j,i),m,n,np,nm) |> transpose
    end
    A
end

function get_lower_right_mat(F,n,m,nm,np)
    getmat(F,2n,2m,nm,np)[end-n+1:end,end-m+1:end]
end

function get_upper_left_mat(F,n,m,nm,np)
    getmat(F,2n,2m,nm,np)[1:n,1:m]
end

function get_upper_mat(F,n,m,nm,np)
    getmat(F,2n,m,nm,np)[1:n,:]
end

function get_lower_mat(F,n,m,nm,np)
    getmat(F,2n,m,nm,np)[n+1:end,:]
end

## ℕ₊ → ℕ₊
function Matrix(Op::SemiLazyBandedOperator{T, S},n,m) where  {T <: ℕ₊, S <: ℕ₊}
    if n > size(Op.A)[1] || m > size(Op.A)[2]
        Op.A = Op.mat(max(n,m))
    end
    Op.A[1:n,1:m]
end

function Matrix(Op::BasicBandedOperator{T, S},n,m) where  {T <: ℕ₊, S <: ℕ₊}
    # or getlrmat(Op.A,n,m,Op.nm,Op.np)
    A = spzeros(n,m)
    if Op.nm > n
        nm = n
    else
        nm = Op.nm
    end
    
    if Op.np > m
        np = m
    else
        np = Op.np
    end
    
    for i = -nm:-1
        if -i <= n-1
            A += spdiagm(n,m, (i => [Op.A(1 -i + j,1 + j) for j in 0:min(n+i-1,m-1) ]))
        end
    end
    
    for i = 1:np
        if i <= m-1
            A += spdiagm(n,m, (i => [Op.A(1 + j,1 + j + i) for j in 0:min(m-i-1,n-1)]))
        end
    end
    
    if Op.nm >= 0 && Op.np >= 0
        A += spdiagm(n,m, [Op.A(i,i) for i in 1:min(n,m)] )
    end
    A
end

## ℕ₋ → ℕ₋
function Matrix(Op::SemiLazyBandedOperator{T, S},n,m) where  {T <: ℕ₋, S <: ℕ₋}
    if n > size(Op.A)[1] || m > size(Op.A)[2]
        Op.A = Op.mat(max(n,m))
    end
    Op.A[end-n+1:end,end-m+1:end]
end

function Matrix(Op::BasicBandedOperator{T, S},n,m) where  {T <: ℕ₋, S <: ℕ₋}  ## Not yet tested
    get_upper_left_mat(Op.A,n,m,Op.nm,Op.np)
end

## ℤ → ℤ
function Matrix(Op::SemiLazyBandedOperator{T,S},n,m) where  {T <: ℤ, S <: ℤ}
    if n > size(Op.A)[1] || m > size(Op.A)[2]
        Op.A = Op.mat(max(n,m))
    end
    N,M = size(Op.A)
    Mm = N₋(M)
    mm = N₋(m)
    Nm = N₋(N)
    nm = N₋(n)
    Op.A[Nm - nm + 1: Nm - nm + 1 + n, Mm - mm + 1: Mm - mm + 1 + m]
end

function Matrix(Op::BasicBandedOperator{T,S},n,m) where {T <: ℤ, S <: ℤ}
    getmat(Op.A,n,m,Op.nm,Op.np)
end

## ℤ → ℕ₋
function Matrix(Op::BasicBandedOperator{T,S},n,m) where {T <: ℤ, S <: ℕ₋}
    get_upper_mat(Op.A,n,m,Op.nm,Op.np)
end

## ℤ → ℕ₊
function Matrix(Op::BasicBandedOperator{T,S},n,m) where {T <: ℤ, S <: ℕ₊}
    get_lower_mat(Op.A,n,m,Op.nm,Op.np)
end

## TODO: Other cases as needed

function Matrix(Op::ProductOfBandedOperators{T,S},n,m) where  {T, S}
    cols = m
    rows = max(cols + rowgrowth(Op.V[end]),1)
    A = Matrix(Op.V[end],rows,cols)
    for j = length(Op.V)-1:-1:2
        cols = rows
        rows = max(cols + rowgrowth(Op.V[j]),1)
        A = Matrix(Op.V[j],rows,cols)*A
    end
    cols = rows
    rows = n
    Matrix(Op.V[1],n,cols)*A
end

function *(Op::ProductOfBandedOperators{T,S},c::Vector) where {T, S}
    cols = length(c)
    rows = max(cols + rowgrowth(Op.V[end]),1)
    v = Matrix(Op.V[end],rows,cols)*c
    for j = length(Op.V)-1:-1:1
        cols = rows
        rows = max(cols + rowgrowth(Op.V[j]),1)
        v = Matrix(Op.V[j],rows,cols)*v
    end
    v
end


