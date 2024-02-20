struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector
end

### NEEDS TO BE FINISHED ###
struct Transform{T<:Basis}
    basis::T
end

function (tf::Transform)(f::BasisExpansion)
    if isconvertible(f.basis,T.basis)
        return Conversion(T.basis)*f
    else
        @error "Not convertible.  Need to discretize."
    end
end

struct plan_transform{T<:Basis}
    basis::T
    n::Integer # redundant
    grid::Vector
    tocoef::Function
end

function *(tf::plan_transform,f::Function)
    tf.tocoef(f.(tf.grid))
end

function *(tf::plan_transform,v::Vector)
    tf.tocoef(v)
end

function *(tf::plan_transform,f::BasisExpansion)
    if isconvertible(f.basis,tf.basis)
       g = Conversion(tf.basis)*f
       return BasisExpansion(g.basis,g[1:tf.n])
    else
        return tf*(x -> f(x))
    end
end

function *(tf::Transform,v::Vector)
    tf(length(v))*v
end

# function BasisExpansion(f::Function,b::Basis,n::Integer)
#     tf = get_transform(b.Domain,b)
#     BasisExpansion(b,tf(n)*f)
# end
####################

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