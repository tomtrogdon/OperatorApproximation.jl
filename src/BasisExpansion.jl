struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector
end

function plot(f::BasisExpansion;dx = 0.01)
    x = -1:dx:1
    x = f.basis.GD.D.map.(x)
    plot(x,f.(x))
end

function plot!(f::BasisExpansion;dx = 0.01)
    x = -1:dx:1
    x = f.basis.GD.D.map.(x)
    plot!(x,f.(x))
end

function *(Op::AbstractOperator,f::BasisExpansion)
    Opc = Op*f.basis
    Opc*f
end

function *(Op::ConcreteLazyOperator,f::BasisExpansion)
    BasisExpansion(Op.range,Op.L*f.c)
end

function pad(v,n)
    if length(v) == n
        return v
    elseif length(v) < n
        return vcat(v,zeros(typeof(v[1]),n-length(v)))
    else
        return v[1:n]
    end
end

function +(f::BasisExpansion,g::BasisExpansion)
    n = max(length(f.c), length(g.c))
    BasisExpansion(f.basis, pad(f.c,n) + pad(g.c,n))
end

function norm(f::BasisExpansion)
    norm(f.c)
end

function *(c::Number,f::BasisExpansion)
    BasisExpansion(f.basis, c*f.c)
end