domainplot(f::BasisExpansion{T};kwargs...) where T = domainplot(f.basis.bases;kwargs...)
coefplot(f::BasisExpansion{T};kwargs...) where T = plot(abs.(f.c) .+ eps(); yaxis = :log, kwargs...)
coefplot!(f::BasisExpansion{T};kwargs...) where T = plot!(abs.(f.c) .+ eps(); yaxis = :log, kwargs...)
coefplot!(p,f::BasisExpansion{T};kwargs...) where T = plot!(p,abs.(f.c) .+ eps(); yaxis = :log, kwargs...)

domainplot(b::DirectSum) = domainplot(b.bases)

domainplot(rhp::RHP;kwargs...) = domainplot(rhdomain(rhp.Γ);kwargs...)

domainplot_kwargs = (aspectratio = 1, arrow = true, linecolor = :dodgerblue, width = 5, legend = false)

function domainplot(D::Interval;kwargs...)
    xs = -1.0:0.001:1.0
    zs = D.map.(xs)
    endpts = [zs[1],zs[end]]
    plot(real(zs),imag(zs); domainplot_kwargs...,kwargs...)
    scatter!(real(endpts),imag(endpts), markersize = 5, markercolor = :black; kwargs...)
end

function domainplot(D::DiscreteDomain;kwargs...)
    scatter(real(D.pts),imag(D.pts), markersize = 5, markercolor = :black; kwargs...)
end

function domainplot!(D::DiscreteDomain;kwargs...)
    scatter!(real(D.pts),imag(D.pts), markersize = 5, markercolor = :black; kwargs...)
end

function domainplot!(D::Interval;kwargs...)
    xs = -1.0:0.001:1.0
    zs = D.map.(xs)
    endpts = [zs[1],zs[end]]
    plot!(real(zs),imag(zs); domainplot_kwargs...,kwargs...)
    scatter!(real(endpts),imag(endpts), markersize = 5, markercolor = :black; kwargs...)
end

function domainplot!(p,D::Interval;kwargs...)
    xs = -1.0:0.001:1.0
    zs = D.map.(xs)
    endpts = [zs[1],zs[end]]
    plot!(p,real(zs),imag(zs); domainplot_kwargs...,kwargs...)
    scatter!(p,real(endpts),imag(endpts), markersize = 5, markercolor = :black; kwargs...)
end

function domainplot(D::Circle;kwargs...)
    zs = (-1.0:0.001:1.0)*pi
    zs = exp.(1im*zs)
    zs = D.map.(zs)
    endpts = [zs[1],zs[end]]
    plot(real(zs),imag(zs); domainplot_kwargs...,kwargs...)
    scatter!(real(endpts),imag(endpts), markersize = 5, markercolor = :black; kwargs...)
end

function domainplot!(D::Circle;kwargs...)
    zs = (-1.0:0.001:1.0)*pi
    zs = exp.(1im*zs)
    zs = D.map.(zs)
    endpts = [zs[1],zs[end]]
    plot!(real(zs),imag(zs); domainplot_kwargs...,kwargs...)
    scatter!(real(endpts),imag(endpts), markersize = 5, markercolor = :black; kwargs...)
end

function domainplot!(p,D::Circle;kwargs...)
    zs = (-1.0:0.001:1.0)*pi
    zs = exp.(1im*zs)
    zs = D.map.(zs)
    endpts = [zs[1],zs[end]]
    plot!(p,real(zs),imag(zs); domainplot_kwargs...,kwargs...)
    scatter!(p,real(endpts),imag(endpts), markersize = 5, markercolor = :black; kwargs...)
end

domainplot(GD::GridDomain;kwargs...) = domainplot(GD.D;kwargs...)
domainplot!(GD::GridDomain;kwargs...) = domainplot!(GD.D;kwargs...)
domainplot(V::Vector{T};kwargs...) where T <: GridDomain = domainplot([GD.D for GD in V];kwargs...)

domainplot(b::Basis;kwargs...) = domainplot(b.GD.D;kwargs...)
domainplot!(b::Basis;kwargs...) = domainplot!(b.GD.D;kwargs...)
domainplot(V::Vector{T};kwargs...) where T <: Basis = domainplot([b.GD.D for b in V];kwargs...)

function domainplot(v::Vector{T};kwargs...) where T <: Domain
    p = domainplot(v[1];kwargs...)
    for i = 2:length(v)
        domainplot!(p, v[i];kwargs...)
    end
    p
end

function coefplot(f::BasisExpansion{T};kwargs...) where T <: DirectSum
    p = coefplot(f[1];kwargs...)
    for i = 2:length(f)
        coefplot!(p,f[i];kwargs...)
    end
    p
end

### needs to be extended
function plot(f::BasisExpansion{T};dx = 0.01,L = 10,kwargs...) where T <: Union{Hermite,OscRational}
    x = -L:dx:L
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    a = f.basis.GD.D.map.(-L)
    b = f.basis.GD.D.map.(L)
    if isreal(a) && isreal(b)# && a < b
        plot(x, y |> real;kwargs...)
        plot!(x, y |> imag;kwargs...)
    else # plot according to arclength
        plot(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function plot!(f::BasisExpansion{T};dx = 0.01,L = 10,kwargs...) where T <: Union{Hermite,OscRational}
    x = -L:dx:L
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    a = f.basis.GD.D.map.(-L)
    b = f.basis.GD.D.map.(L)
    if isreal(a) && isreal(b)# && a < b
        plot!(x, y |> real;kwargs...)
        plot!(x, y |> imag;kwargs...)
    else # plot according to arclength
        plot!(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function plot(f::BasisExpansion{T};dx = 0.01,L = 10,kwargs...) where T <: Laguerre
    x = 0:dx:L
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    a = f.basis.GD.D.map.(0)
    b = f.basis.GD.D.map.(L)
    display(x)
    if isreal(a) && isreal(b)# && a < b
        plot(x |> real, y |> real;kwargs...)
        plot!(x |> real, y |> imag;kwargs...)
    else # plot according to arclength
        plot(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function plot!(f::BasisExpansion{T};dx = 0.01,L = 10,kwargs...) where T <: Laguerre
    x = 0:dx:L
    x = f.basis.GD.D.map.(x)
    y = f.(x)
    a = f.basis.GD.D.map.(0)
    b = f.basis.GD.D.map.(L)
    if isreal(a) && isreal(b)# && a < b
        plot!(x, y |> real;kwargs...)
        plot!(x, y |> imag;kwargs...)
    else # plot according to arclength
        plot!(abs.(x .- a), y |> real;kwargs...)
        plot!(abs.(x .- a), y |> imag;kwargs...)
    end
end

function plot(f::BasisExpansion{T};dx = 0.01,kwargs...) where T
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

function plot!(f::BasisExpansion{T};dx = 0.01,kwargs...) where T
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

function rhplot(rhp::RHP;kwargs...)
    # need to extend for larger RHPs
    dom = rhp.Γ |> rhdomain
    p0 = domainplot(dom;kwargs...)
    if length(rhp.P) > 0
        scatter!(p0,real(rhp.P),imag(rhp.P);markercolor = :lightblue,kwargs...)
    end
    ran = rhrange(dom)
    N = 100
    plts = [p0]
    y = 0
    for i = 1:size(rhp.Γ,1)
        d = dom[i]
        x = d.GD.grid(N)
        z = d.GD.D.map.(x)
        y = vcat(map( x -> reshape(mofeval(rhp.J[i],x),1,:), z)...)
        p1 = plot(x,y[:,1] |> real;legend = false, kwargs...)
        plot!(p1,x,y[:,1] |> imag;legend = false, kwargs...)
        for k = 2:size(y,2)
            plot!(p1,x,y[:,k] |> real;legend = false, kwargs...)
            plot!(p1,x,y[:,k] |> imag;legend = false, kwargs...)
        end
        push!(plts,p1)
    end
    plts
end

function ploteval(E::ContinuousEigen)
    scatter(E.values |> real, E.values |> imag)
end

function ploteval!(E::ContinuousEigen)
    scatter!(E.values |> real, E.values |> imag)
end