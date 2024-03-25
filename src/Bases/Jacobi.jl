struct Jacobi <: Basis
    α::Number
    β::Number
    GD::GridInterval
end
####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::Jacobi) = ℕ₊

function dim(sp::Jacobi)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Jacobi
    testconv(f.c)
end

function chop(f::BasisExpansion{T}) where T <: Jacobi
    BasisExpansion(f.basis,chop(f.c))
end
####################################
####################################
####################################

Tgrid = n -> cos.( (2*(1:n) .- 1)/(2*n) * pi ) |> reverse

function Jacobi_ab(a,b) #TODO: simplify evaluation
    bfun = n -> (a+b==-1 && n==0) ? √(2*a*b) :
        2*sqrt(n+1)*sqrt(n+a+1)*sqrt(n+b+1)*sqrt(n+a+b+1)/
        ((2n + a +b + 2)*sqrt(2n + a +b + 3)*sqrt(2n + a +b + 1))
    afun = n -> ((a+b==0 || a+b==-1) && n==0) ? (b-a)/(a+b+2) :
        (b^2 - a^2)/((2n + a +b + 2)*(2n + a +b))
    return (n -> n < 0 ? a : afun(n),n -> n < 0 ? b : bfun(n))  # this is not needed but lets us get access
                                                                # to the Jacobi parameters after the fact
end

function jacobi(a,b,n) # creates (n + 1) x (n+1) Jacobi matrix
   SymTridiagonal([a(i) for i in 0:n],[b(i) for i in 0:n-1])
end

function Interp_transform(a,b,n)
    E = jacobi(a,b,n) |> eigen
    return E.values, E.vectors*(Diagonal(E.vectors[1,:]))
end

function OPMultiplication(a::Function,b::Function,α::Function,β::Function,c::Vector,u::Union{Array,SparseMatrixCSC})
    ## α, β give the recurrence coefficients for the basis used
    # on the domain of the operator.  a, b give the recurrence coefficitions
    # for the basis used for the function doing the multiplication.       
    n = size(u)[1] + length(c) + 3
    J = jacobi(α,β,n-1)
    #display(J)
    shape = size(u)
    shape = tuple(length(c)+ 3,shape[2:end]...)
    U = vcat(u,zeros(shape)) |> sparse # TODO: Better way to do this?
    p = U # p_0
    q = c[1]*p
    polder = p
    p = J*p - a(0)*p # compute p_1
    p /= b(0)
    q += c[2]*p
    pold = p
    for j = 3:length(c) # compute p_n
        p = J*pold - a(j-2)*pold - b(j-3)*polder
        p /= b(j-2)
        q += c[j]*p
        polder = pold
        pold = p
    end
    q |> sparse
end

function Gauss_quad(a,b,n)
    E = jacobi(a,b,n) |> eigen
    return E.values, abs2.(E.vectors[1,:])
end

function Interp_transform(f::Function,a,b,n)
    E = jacobi(a,b,n) |> eigen
    E.vectors*(Diagonal(E.vectors[1,:])*map(f,E.values))
end

function poly(a,b,n,x) # a naive use of the three-term recurrence
    if n == 1
        return fill(1.0,n)
    end
    p = fill(0.0,n)
    p[1] = 1.0 # p_0
    p[2] = x.*p[1] - a(0)*p[1] # compute p_1
    p[2] /= b(0)
    for j = 1:n-2 # compute p_n
        p[j+2] = x.*p[j+1] - a(j)*p[j+1] - b(j-1)*p[j]
        p[j+2] /= b(j)
    end
    p
end

function poly(a,b,n,z::Vector)
    vcat(map(zz -> poly(a,b,n,zz) |> transpose , z)...)
end

function JacobiW(a,b,x)
    c = (_₂F₁(1, -a, 2+b, -1))/(1+b) + (_₂F₁(1, -b, 2+a, -1))/(1+a)
    return 1/c*(1-x)^a*(1+x)^b
end

function JacobiWconst(a,b)
    c = (_₂F₁(1, -a, 2+b, -1))/(1+b) + (_₂F₁(1, -b, 2+a, -1))/(1+a)
    return 1/c
end

normalization(n::Int,α::Real,β::Real) = 2^(α+β)*gamma(n+α+1)*gamma(n+β+1)/gamma(2n+α+β+2)
#stieltjesjacobimoment(α::Real,β::Real,n::Int,z) =
    #(x = 2/(1-z);normalization(n,α,β)*HypergeometricFunctions.mxa_₂F₁(n+1,n+α+1,2n+α+β+2,x))
stieltjesjacobimoment(α::Real,β::Real,n::Int,z) =
    (x = 2/(1-z);HypergeometricFunctions.mxa_₂F₁(n+1,n+α+1,2n+α+β+2,x))/2
stieltjesjacobimoment(α::Real,β::Real,z) = stieltjesjacobimoment(α,β,0,z)
JacobiSeed(α,β,z) = 1im/(2*pi)*stieltjesjacobimoment(α,β,z)