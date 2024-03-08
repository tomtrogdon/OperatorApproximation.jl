struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector # if DirectSum then c is a vector of vectors
end

function BasisExpansion(f::Function,basis::Basis,N::Integer)
    Conversion(basis)*BasisExpansion(f,GridValues(basis.GD),N)
end

function BasisExpansion(f::BasisExpansion,sp::Basis)
    Conversion(sp)*f
end

# This is, in effect, a projection.
function BasisExpansion(f::BasisExpansion,sp::Basis,N::Integer)
    # Would produce incorrect results for DiscreteBasis input
    if typeof(f.basis) <: DiscreteBasis
        @error "Unable to project a DiscreteBasis"
        return
    end
    display(sp)
    g = Conversion(sp)*f
    if length(g.c) < N  # TODO:  Not correct for bi-infinite vectors
        @warn "Input dimension smaller than linear system size. Padding with zeros."
    end
    BasisExpansion(g.basis,pad(g,N))
end

### needs to be extended
function plot(f::BasisExpansion;dx = 0.01)
    x = -1:dx:1
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    plot(x, y |> real)
    plot!(x, y |> imag)
end

function plot!(f::BasisExpansion;dx = 0.01)
    x = -1:dx:1
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    plot!(x, y |> real)
    plot!(x, y |> imag)
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

function chop(c::Vector)
    ind = length(c)
    for i = length(c):-1:1
        if norm(c[i:end]) > 1e-15
            ind = i
            break
        end
    end
    c[1:ind]
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