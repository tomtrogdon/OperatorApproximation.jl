struct MarchenkoPastur <: Basis
    d::Number
    GD::GridInterval
end

####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::MarchenkoPastur) = ℕ₊

function dim(sp::MarchenkoPastur)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: MarchenkoPastur
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: MarchenkoPastur
    BasisExpansion(f.basis,chop(f.c))
end

function getweight(sp::MarchenkoPastur)
    d = sp.d
    x -> 2/(pi*(2*sqrt(d)*x + 1 + d))*sqrt(1-x)*sqrt(1 + x)
end
####################################
#####  Important to implement  #####
####################################
function sum(f::BasisExpansion{T}) where T <: MarchenkoPastur
    (f.basis.GD.D.b - f.basis.GD.D.a)*f.c[1]
end

function moment(f::BasisExpansion{T},k::Int64) where T <: MarchenkoPastur
    if k == 0
        return sum(f)
    end
    Multiplication(x -> x^k)*f |> sum
end

function (P::BasisExpansion{MarchenkoPastur})(X::Number) # Clenshaw's algorithm
    n = P.c |> length
    d = P.basis.d
    x = P.basis.GD.D.imap(X)
    a,b = MP_ab(d)
    (hcat(e(1,n) |> sparse,(jacobi(a,b,n) - x*I)[1:end-1,1:end-2] |> sparse)\P.c)[1]
end