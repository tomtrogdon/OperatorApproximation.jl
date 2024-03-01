struct Ultraspherical <: Basis
    λ::Number
    GD::GridDomain
end

function Chop(f::BasisExpansion{T}) where T <: Ultraspherical
    ind = length(f.c)
    for i = length(f.c):-1:1
        if norm(f.c[i:end]) > 1e-15
            ind = i
            break
        end
    end
    BasisExpansion(f.basis,f.c[1:ind])
end

function dim(sp::Ultraspherical)
    Inf
end

function (P::BasisExpansion{Ultraspherical})(X::Number) # Clenshaw's algorithm
    n = P.c |> length
    λ = P.basis.λ
    x = P.basis.GD.D.imap(X)
    a,b = Jacobi_ab(λ - 1/2, λ - 1/2)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
end

# function BasisExpansion(f::Function,basis::Ultraspherical,N::Integer)
#     Conversion(basis)*BasisExpansion(f,GridValues(basis.GD),N)
# end
