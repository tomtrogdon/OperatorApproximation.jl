struct Hardy{T <: GridRegion, S <: Domain} <: Basis
    GD::T
end
Hardy(GD) = Hardy{typeof(GD),typeof(GD.GD.D)}(GD)

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::Hardy{T,S})  where {T <: Exterior, S <: Circle} = ℕ₋
cfd(sp::Hardy{T,S})  where {T <: Interior, S <: Circle} = ℕ₊

function dim(sp::Hardy)
    Inf
end

function pad(f::BasisExpansion{T},N) where T <: Hardy
    BasisExpansion(f.basis,pad(f.c,N))
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

function (P::BasisExpansion{Hardy{T,S}})(X::Number) where {T <: Exterior, S <: Circle}
    n = P.c |> length
    x = 1.0./P.basis.GD.D.imap(X)
    if abs(x) > 1
        return 0.0im
    end
    p = x
    sum = P.c[1]*p
    for i = 2:n
        p *= x
        sum += P.c[i]*p        
    end
    sum
end

function (P::BasisExpansion{Hardy{T,S}})(X::Number) where {T <: Interior, S <: GridCircle}
    n = P.c |> length
    x = P.basis.GD.D.imap(X)
    if abs(x) > 1
        return 0.0im
    end
    p = 1.0
    sum = P.c[1]
    for i = 2:n
        p *= x
        sum += P.c[i]*p        
    end
    sum
end