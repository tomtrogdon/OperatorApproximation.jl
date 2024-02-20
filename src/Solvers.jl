function \(L::Vector{ConcreteOperator},b::Vector,N::Integer)
    k = 0
    Ops = []
    rhss = []
    for i = 1:length(L)
        if isfinite(rank(L[i]))
            k += rank(L[i])
            push!(Ops,Matrix(L[i],N))
            push!(rhss,b[i])
        end
    end
    Op = vcat(Ops...)
    rhs = vcat(rhss...)
    ## Assume only one non Functional Operator, for now...
    ## Adjust # of rows if not...
    for i = 1:length(L)
        if !isfinite(rank(L[i]))
            Op = vcat(Op,Matrix(L[i],N-k,N))
            T = Transform(L[i].range)
            rhs = vcat(rhs,T(N-k)*b[i])
        end
    end
    BasisExpansion(L[1].domain,Op\rhs)
end

function testconv(f::BasisExpansion)
    norm(f.c[end-4:end]) < tol
end


function \(L::Vector{ConcreteOperator},b::Vector)
    if !(typeof(N) <: Integer)
        n = 32
        sol = \(L,b,n)
        bool = testconv(sol)
        while !bool
            n *= 2
            sol = \(L,b,n)
            bool = testconv(sol)
        end
        return sol
    else
        return \(L,b,N)
    end
end

function \(L::Vector{AbstractOperator},b::Vector)
    Ops = [J*basis for J in L]
    Ops\b
end

function \(L::Vector{AbstractOperator},b::Vector,basis::Basis,N::Integer)
    Ops = [J*basis for J in L]
    \(Ops,b,N)
end

function \(L::Vector{AbstractOperator},b::Vector,N::Integer)
    Ops = [J*basis for J in L]
    \(Ops,b,N)
end