struct Hardy{T <: GridRegion, S <: Domain} <: Basis
    GD::T
end
Hardy(GD) = Hardy{typeof(GD),typeof(GD.GD.D)}(GD)

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::Hardy{T,S})  where {T <: Exterior, S <: Circle} = ℕ₋
cfd(sp::Hardy{T,S})  where {T <: Interior, S <: Circle} = ℕ₊
cfd(sp::Hardy{T,S})  where {T <: Exterior, S <: DiscreteDomain} = ℕ₊
cfd(sp::Hardy{T,S})  where {T <: Exterior, S <: Interval} = ℕ₊
cfd(sp::Hardy{T,S})  where {T <: Exterior, S <: Axis} = ℕ₋
cfd(sp::Hardy{T,S})  where {T <: Interior, S <: Axis} = ℕ₊



function dim(sp::Hardy)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Hardy
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: Hardy  # add tolerance?
    BasisExpansion(f.basis,chop(f.c))
end
####################################
####################################
####################################

function p_partial_horner(c,x)
    n = length(c)
    if abs(x) > 1 || isnan(x)
        return 0.0im
    end
    p = 1.0
    sum = c[1]
    for i = 2:n
        p *= x
        sum += c[i]*p        
    end
    sum
end

function n_partial_horner(c,x)
    n = length(c)
    if abs(x) >= 1 || isnan(x)
        return 0.0im
    end
    p = x
    sum = c[end]*p
    for i = n-1:-1:1
        p *= x
        sum += c[i]*p        
    end
    sum
end

function (P::BasisExpansion{Hardy{T,S}})(X::Number) where {T <: Exterior, S <: Circle}
    n_partial_horner(P.c,1.0./P.basis.GD.D.imap(X))
end

function (P::BasisExpansion{Hardy{T,S}})(X::Number) where {T <: Interior, S <: Circle}
    p_partial_horner(P.c,P.basis.GD.D.imap(X))
end

function (P::BasisExpansion{Hardy{T,S}})(X::Number) where {T <: Interior, S <: Axis}
    x = P.basis.GD.D.imap(X)
    if imag(x) < 0
        return 0.0im
    end
    x = ((x.-1im)./(x.+1im))
    p_partial_horner(P.c,x)*x - sum(P.c)*(abs(x) <= 1 ? 1.0 : 0)
end

function (P::BasisExpansion{Hardy{T,S}})(X::Number) where {T <: Exterior, S <: Axis}
    x = P.basis.GD.D.imap(X)
    if imag(x) >= 0
        return 0.0im
    end
    x = ((x.+1im)./(x.-1im))
    n_partial_horner(P.c,x) - sum(P.c)*(abs(x) >= 1 ? 0.0 : 1.0)
end

function (P::BasisExpansion{Hardy{Exterior{T},S}})(X::Number) where {T <: Union{JacobiMappedInterval,JacobiInterval}, S <: Interval}
    α = P.basis.GD.GD.α
    β = P.basis.GD.GD.β
    a, b = Jacobi_ab(α,β)
    dot(cauchy(a,b,JacobiSeed(α,β),length(P.c)-1,P.basis.GD.GD.D.imap(X)) |> conj,P.c)*2
end

function (P::BasisExpansion{Hardy{Exterior{T},S}})(X::Number) where {T <: Union{MarchenkoPasturMappedInterval,MarchenkoPasturInterval}, S <: Interval}
    d = P.basis.GD.GD.d
    a, b = MP_ab(d)
    dot(cauchy(a,b,MPSeed(d),length(P.c)-1,P.basis.GD.GD.D.imap(X)) |> conj,P.c)*2
end

function (P::BasisExpansion{Hardy{T,S}})(X::Number) where {T <: Exterior, S <: DiscreteDomain}
    dot(poleres_cauchy(P.basis.GD.grid,X) |> conj,P.c)
end