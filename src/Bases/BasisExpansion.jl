struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector # if DirectSum then c is a vector of vectors
end

isDirectSum(f::BasisExpansion) = typeof(f.c[1]) <: Vector

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

function Base.chop(c::Vector)
    ind = length(c)
    for i = length(c):-1:1
        if norm(c[i:end]) > 1e-15
            ind = i
            break
        end
    end
    c[1:ind]
end

function Base.chop(f::BasisExpansion{T}) where T
    if cfd(f.basis) == ℤ
        nm = N₋(length(f.c))
        vm = Base.chop(copy(f.c[1:nm]) |> reverse) |> reverse;
        vp = Base.chop(copy(f.c[nm+1:end]));
        k = max(length(vm),length(vp))
        vm = pad(vm |> reverse,k) |> reverse
        vp = pad(vp,k)
        return BasisExpansion(f.basis,vcat(vm,vp))
    elseif cfd(f.basis) == ℕ₊
        coeffs = Base.chop(copy(f.c));
        return BasisExpansion(f.basis,coeffs)
    elseif cdf(f.basis) == ℕ₋
        coeffs = Base.chop(copy(f.c) |> reverse) |> reverse;
        return BasisExpansion(f.basis,coeffs)
    else
        @warn "chop() hasn't been implemented for this coefficient domain yet"
    end
end

function Base.chop(f::BasisExpansion{T}) where T <: DirectSum
    dummy = Base.chop(f[1]) ⊕ Base.chop(f[2])
    if length(f) == 2
        return dummy 
    else
        for i=3:length(f)
            dummy = dummy ⊕ Base.chop(f[i])
        end
        return dummy
    end
end

function Base.chop(v::Vector{T}) where T <: BasisExpansion
    for i=1:length(v)
        v[i] = Base.chop(v[i])
    end
    return v
end

function combinebasexp(f::Vector{T}) where T <: BasisExpansion
    j = 0
    @inbounds while j < length(f)
        j += 1
        if maximum(abs.(f[j].c)) < 1e-15
            deleteat!(f,j)
            j -= 1
        end
    end
    i = 0
    @inbounds while i < length(f)
        i += 1
        j = i
        while j < length(f)
            j += 1
            if isconvertible(f[i].basis,f[j].basis)
                f[i] += f[j]
                deleteat!(f,j)
                j -= 1
            end
        end
    end
    if length(f) == 1
        return f[1]
    else
        return f
    end
end

function combine(f::Vector{T}) where T <: BasisExpansion
    count = 0
    for i=1:length(f)
        if typeof(f[i]) == BasisExpansion{DirectSum}
            f[i] = combine(f[i])
            if typeof(f[i]) == BasisExpansion{DirectSum}
                count += 1
            end
        end
    end
    testvec = Vector{Any}(undef,count)
    testvec2 = Vector{BasisExpansion}(undef,length(f)-count)
    ind = 1
    ind2 = 1
    for i=1:length(f)
        if typeof(f[i]) == BasisExpansion{DirectSum}
            testvec[ind] = f[i]
        else
            testvec2[ind2] = f[i]
        end
    end
    BasExp = combinebasexp(testvec2)
    if isempty(testvec)
        return BasExp
    elseif typeof(BasExp) == BasisExpansion
        return push!(testvec,BasExp)
    else
        return vcat(testvec,BasExp)
    end
end

function combine(f::BasisExpansion{T}) where T <: DirectSum
    j = 0
    @inbounds while j < length(f)
        j += 1
        if maximum(abs.(f[j].c)) < 1e-15
            if j == 1
                f = f[2:end]
            elseif j == length(f)
                f = f[1:end-1]
            elseif j == 2
                dummy = f[1] ⊕ f[3]
                for k=4:length(f)
                    dummy = dummy ⊕ f[k]
                end
                f = dummy
            else
                dummy = f[1] ⊕ f[2]
                for k=3:length(f)
                    if k==j
                        continue
                    else
                        dummy = dummy ⊕ f[k]
                    end
                end
                f = dummy  
            end
            j -= 1
        end
    end

    array = collect(1:length(f))
    testvec = Vector{Any}(undef,length(f))
    i = 0
    count = 0
    @inbounds while i < length(array)
        i += 1
        count += 1
        j = i 
        testvec[count] = f[array[i]]
        while j < length(array)
            j += 1
            if isconvertible(f[array[i]].basis,f[array[j]].basis)
                testvec[count] = testvec[count] + f[array[j]]
                deleteat!(array,j)
                j -= 1
            end
        end
    end
    dsums = testvec[1:count]
    if length(dsums) == 1
        return dsums[1]
    else
        g = dsums[1] ⊕ dsums[2]
        for i=3:length(dsums)
            g = g ⊕ dsums[i]
        end
        return g
    end
end

function combine(f::Vector{BasisExpansion{T}}) where T <: DirectSum
    testvec = Vector{BasisExpansion}(undef,length(f))
    for i=1:length(f)
        testvec[i] = combine(f[i])
    end
    return testvec
end

function simp(f::Vector{T}) where T <: BasisExpansion
    return chop(combine(chop(f)))
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
    if length(f.basis) == length(g.basis) && f.basis == g.basis
        out = []
        for i = 1:length(f.c)
            n = max(length(f.c[i]), length(g.c[i]))
            tp = cfd(f.basis.bases[i])
            push!(out, pad(tp,f.c[i],n) + pad(tp,g.c[i],n))
        end
        BasisExpansion(f.basis,out)
    else
        @warn "Basis not compatible, using direct sum"
        f ⊕ g
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