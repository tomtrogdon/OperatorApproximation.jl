struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector
end

function BasisExpansion(f::BasisExpansion,sp::Basis)
    Conversion(sp)*f
end

# This is, in effect, a projection.
function BasisExpansion(f::BasisExpansion,sp::Basis,N::Integer)
    # Would produce incorrect results for DiscreteBasis input
    if f.basis <: DiscreteBasis
        @error "Unable to project a DiscreteBasis"
        return
    end
    g = Conversion(sp)*f
    if g.c < N
        @error "Input dimension too small"
        return
    end
    BasisExpansion(g.basis,g.c[1:N])
end

### needs to be extended

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