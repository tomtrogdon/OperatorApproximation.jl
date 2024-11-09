struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector # if DirectSum then c is a vector of vectors
end

transpose(f::BasisExpansion{T}) where T = f

function (V::Vector{T})(z) where T <: BasisExpansion
    [v(z) for v in V]
end

function (V::Matrix{T})(z) where T <: BasisExpansion
    [v(z) for v in V]
end

function sum(f::BasisExpansion{T}) where T<:DirectSum
    sum([sum(f[i]) for i = 1:length(f)])
end

function moment(f::BasisExpansion{T},k) where T<:DirectSum
    sum([moment(f[i],k) for i = 1:length(f)])
end

getindex(b::BasisExpansion{T},i::Int64) where T <: DirectSum = BasisExpansion(b.basis[i],b.c[i])
getindex(b::BasisExpansion{T},i::UnitRange{Int64}) where T <: DirectSum = BasisExpansion(b.basis[i],b.c[i])
getindex(b::BasisExpansion{T},i) where T = getindex([b],i)
getindex(b::BasisExpansion{T},i::UnitRange{Int64}) where T = i == 1:1 ? b : @error "out of range"
axes(b::BasisExpansion{T}) where T <: DirectSum = axes(b.basis.bases)
axes(b::BasisExpansion{T},i) where T <: DirectSum = axes(b.basis.bases,i)
axes(b::BasisExpansion{T}) where T = axes([b])
axes(b::BasisExpansion{T},i) where T = axes([b],i)
size(b::BasisExpansion{T}) where T <: DirectSum = size(b.basis.bases)
size(b::BasisExpansion{T},i) where T <: DirectSum = size(b.basis.bases,i)
lastindex(b::BasisExpansion{T}) where T <: DirectSum = b.basis.bases |> length
length(b::BasisExpansion{T}) where T <: DirectSum = b.basis.bases |> length
domainplot(f::BasisExpansion{T};kwargs...) where T = domainplot(f.basis.bases;kwargs...)
coefplot(f::BasisExpansion{T};kwargs...) where T = plot(abs.(f.c) .+ eps(); yaxis = :log, kwargs...)
coefplot!(f::BasisExpansion{T};kwargs...) where T = plot!(abs.(f.c) .+ eps(); yaxis = :log, kwargs...)
coefplot!(p,f::BasisExpansion{T};kwargs...) where T = plot!(p,abs.(f.c) .+ eps(); yaxis = :log, kwargs...)

function coefplot(f::BasisExpansion{T};kwargs...) where T <: DirectSum
    p = coefplot(f[1];kwargs...)
    for i = 2:length(f)
        coefplot!(p,f[i];kwargs...)
    end
    p
end

function BasisExpansion(f::Function,basis::Basis,N::Integer)
    Conversion(basis)*BasisExpansion(f,GridValues(basis.GD),N)
end

function ⊕(f::BasisExpansion{T}) where T
   f
end

function ⊕(f::BasisExpansion{T},g::BasisExpansion{S}) where {T, S}
    BasisExpansion(f.basis ⊕ g.basis, [f.c,g.c])
end

function ⊕(f::BasisExpansion{T},g::BasisExpansion{S}) where {T <: DirectSum, S}
    v = copy(f.c)
    push!(v,g.c)
    BasisExpansion(f.basis ⊕ g.basis, v)
end

function ⊕(f::BasisExpansion{T},g::BasisExpansion{S}) where {T, S <: DirectSum}
    v = copy(g.c)
    pushfirst!(v,f.c)
    BasisExpansion(f.basis ⊕ g.basis, v)
end

function ⊕(f::BasisExpansion{T},g::BasisExpansion{S}) where {T <: DirectSum, S <: DirectSum}
    BasisExpansion(f.basis ⊕ g.basis, vcat(f.c,g.c))
end

function BasisExpansion(f::Function,basis::Basis)
    n = 4
    g = Conversion(basis)*BasisExpansion(f,GridValues(basis.GD),n)
    while !testconv(g) && n < Nmax
        n *= 2
        g = Conversion(basis)*BasisExpansion(f,GridValues(basis.GD),n)
    end
    if n >= Nmax  ## make this correct for all bases
        @warn "Max DOF reached.  Tolerance met: "*string(norm(g.c[end-4:end]))
    end
    return g
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
    g = Conversion(sp)*f
    if length(g.c) < N  # TODO:  Not correct for bi-infinite vectors
        @warn "Input dimension smaller than linear system size. Padding with zeros."
    end
    #BasisExpansion(g.basis,pad(g,N).c)
    pad(g,N)
end

### needs to be extended
function plot(f::BasisExpansion;dx = 0.01,kwargs...)
    x = -1:dx:1
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    a = f.basis.GD.D.a
    b = f.basis.GD.D.b
    if isreal(a) && isreal(b) && a < b
        plot(x, y |> real;kwargs...)
        plot!(x, y |> imag;kwargs...)
    else # plot according to arclength
        plot(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function plot!(f::BasisExpansion;dx = 0.01,kwargs...)
    x = -1:dx:1
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    a = f.basis.GD.D.a
    b = f.basis.GD.D.b
    if isreal(a) && isreal(b) && a < b
        plot!(x, y |> real;kwargs...)
        plot!(y, y |> imag;kwargs...)
    else # plot according to arclength
        plot!(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function weightplot(f::BasisExpansion;dx = 0.01,kwargs...)
    X = -1:dx:1
    x = f.basis.GD.D.map.(X)
    y = f.(x)
    w = getweight(f.basis)
    W = w.(X)
    y = W.*y
    a = f.basis.GD.D.a
    b = f.basis.GD.D.b
    y *= 2
    if isreal(a) && isreal(b) && a < b
        plot(x, y |> real;kwargs...)
        plot!(x, y |> imag;kwargs...)
    else # plot according to arclength
        plot(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function weightplot!(f::BasisExpansion;dx = 0.01,kwargs...)
    X = -1:dx:1
    x = f.basis.GD.D.map.(X)
    y = f.(x)
    w = getweight(f.basis)
    W = w.(X)
    y = W.*y
    a = f.basis.GD.D.a
    b = f.basis.GD.D.b
    y *= 2
    if isreal(a) && isreal(b) && a < b
        plot!(x, y |> real;kwargs...)
        plot!(x, y |> imag;kwargs...)
    else # plot according to arclength
        plot!(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function pad(v::Vector,n::Int64)
    if length(v) == n
        return v
    elseif length(v) < n
        return vcat(v,zeros(typeof(v[1]),n-length(v)))
    else
        return v[1:n]
    end
end

function pad(::Type{ℕ₊},v::Vector,n::Int64) 
    pad(v,n)
end

function pad(::Type{ℕ₋},v::Vector,n::Int64)
    pad(v |> reverse ,n) |> reverse
end

function pad(::Type{𝕏},v::Vector,n::Int64)
    v
end

function pad(::Type{ℤ},v::Vector,n::Int64)
    nm = N₋(length(v))
    new_nm = N₋(n)
    vm = pad(ℕ₋,copy(v[1:nm]),new_nm)
    vp = pad(ℕ₊,copy(v[nm+1:end]),n - new_nm)
    vcat(vm,vp)
end

## TODO: Need implmentation for DirectSum
function pad(f::BasisExpansion,N::Int64)
    BasisExpansion(f.basis,pad(cfd(f.basis),f.c,N))
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

function +(f::BasisExpansion{S},g::BasisExpansion{T}) where {S, T}
    if f.basis == g.basis
        tp = cfd(f.basis)
        n = max(length(f.c), length(g.c))
        BasisExpansion(f.basis, pad(tp,f.c,n) + pad(tp,g.c,n))
    else
        @warn "Basis not compatible, using direct sum"
        f ⊕ g
    end
end

function +(f::BasisExpansion{S},g::BasisExpansion{T}) where {S <: DirectSum, T <: DirectSum}
    if f.basis == g.basis
        out = []
        for i = 1:length(f.c)
            n = max(length(f.c[i]), length(g.c[i]))
            tp = cfd(f.basis.bases[i])
            push!(out, pad(tp,f.c[i],n) + pad(tp,g.c[i],n))
        end
        BasisExpansion(f.basis,out)
    else
        @error "Basis not compatible"
    end
end

function norm(f::BasisExpansion{S}) where S
    norm(f.c)
end

function *(c::Number,f::BasisExpansion{S}) where S
    BasisExpansion(f.basis, c*f.c)
end

function -(f::BasisExpansion,g::BasisExpansion)
    f + (-1)*g
end

# Probably inefficient...
function (f::BasisExpansion{T})(x) where T <: DirectSum
    [f[i](x) for i in 1:size(f,1)] |> sum
end