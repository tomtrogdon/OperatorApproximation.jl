struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector
end

struct transform
    grid::Function
    tocoef::Function
end

struct plan_transform
    n::Integer
    grid::Function
    tocoef::Function
end

function (tf::transform)(n)
    plan_transform(n,tf.grid,tf.tocoef)
end

function *(tf::plan_transform,f::Function)
    tf.tocoef(f.(tf.grid(tf.n)))
end

function *(tf::plan_transform,v::vector)
    tf.tocoef(v)
end

function *(tf::transform,v::Vector)
    n = length(v)
    tf(n)*v
end

function BasisExpansion(f::Function,b::Basis,n::Integer)
    tf = get_transform(b.Domain,b)
    BasisExpansion(b,tf(n)*f)
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

function -(f::BasisExpansion,g::BasisExpansion)
    f + (-1)*g
end