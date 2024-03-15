function Matrix(Op::SemiLazyBandedOperator{T, S},n,m) where  {T <: ℕ₊, S <: ℕ₊}
    if n > size(Op.A)[1] || m > size(Op.A)[2]
        Op.A = Op.mat(max(n,m))
    end
    Op.A[1:n,1:m]
end

function Matrix(Op::BasicBandedOperator{T, S},n,m) where  {T <: ℕ₊, S <: ℕ₊}
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

function Matrix(Op::ProductOfBandedOperators{T,S},n,m) where  {T <: ℕ₊, S <: ℕ₊}
    cols = m
    rows = max(cols+Op.V[end].nm,1)
    A = Matrix(Op.V[end],rows,cols)
    for j = length(Op.V)-1:-1:2
        cols = rows
        rows = max(cols + Op.V[j].nm,1)
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
        function *(Op::$op{T,S},c::Vector) where {T <: ℕ₊, S <: ℕ₊}
            n = length(c)
            m = max(Op.nm + n,1)
            Matrix(Op,m,n)*c
        end

        function rowgrowth(Op::$op{T,S}) where {T <: ℕ₊, S <: ℕ₊}
            Op.nm
        end
    end
end

function *(Op::ProductOfBandedOperators{T,S},c::Vector) where {T <: ℕ₊, S <: ℕ₊}
    cols = length(c)
    rows = max(cols+Op.V[end].nm,1)
    v = Matrix(Op.V[end],rows,cols)*c
    for j = length(Op.V)-1:-1:1
        cols = rows
        rows = max(cols + Op.V[j].nm,1)
        v = Matrix(Op.V[j],rows,cols)*v
    end
    v
end

function rowgrowth(Op::ProductOfBandedOperators{T,S}) where {T <: ℕ₊, S <: ℕ₊}
    cols = 0
    rows = cols+Op.V[end].nm
    for j = length(Op.V)-1:-1:1
        cols = rows
        rows = cols + Op.V[j].nm
    end
    rows
end

# function SIIdentityOperator()
#     BasicBandedOperator(SI,0,0, (i,j) -> i == j ? 0.0 : 1)
# end