function Jacobi_ab(a,b) #TODO: simplify evaluation
    bfun = n -> (a+b==-1 && n==0) ? √(2*a*b) :
        2*sqrt(n+1)*sqrt(n+a+1)*sqrt(n+b+1)*sqrt(n+a+b+1)/
        ((2n + a +b + 2)*sqrt(2n + a +b + 3)*sqrt(2n + a +b + 1))
    afun = n -> ((a+b==0 || a+b==-1) && n==0) ? (b-a)/(a+b+2) :
        (b^2 - a^2)/((2n + a +b + 2)*(2n + a +b))
    return (n -> n < 0 ? a : afun(n),n -> n < 0 ? b : bfun(n))  # this is not needed but lets us get access
                                                                # to the Jacobi parameters after the fact
end

function Laguerre_ab(a)
    bfun = n -> sqrt(1 + n)*sqrt(1+ n + a);
    afun = n -> 1 + a + 2n
    return (afun, bfun)
end

function MP_ab(d) # need to map support to [-1,1]
   bfun = n -> 0.5# sqrt(d)
   afun = n -> n == 0 ? -sqrt(d)/2.0 : 0.0
   (afun,bfun)                                   
end

function Hermite_ab()
    bfun = n -> sqrt(n + 1.0)
    afun = n -> 0.0
    (afun, bfun)
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
    shape = size(u)
    shape = tuple(length(c)+ 3,shape[2:end]...)
    U = vcat(u,zeros(shape)) |> sparse # TODO: Better way to do this?
    p = U # p_0
    q = c[1]*p
    if length(c) == 1
        return q
    end
    polder = p
    p = J*p - a(0)*p # compute p_1
    p /= b(0)
    q += c[2]*p
    if length(c) == 2
        return q
    end
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

function MPW(d,x)  # not correct...
    dp = (1 + sqrt(d))^2
    dm = (1 - sqrt(d))^2
    1/(2*pi*d*x)*sqrt(dp-x)*sqrt(x-dm)
end

## from SingularIntegralEquation.jl
normalization(n::Int,α::Real,β::Real) = 2^(α+β)*gamma(n+α+1)*gamma(n+β+1)/gamma(2n+α+β+2)
#stieltjesjacobimoment(α::Real,β::Real,n::Int,z) =
    #(x = 2/(1-z);normalization(n,α,β)*HypergeometricFunctions.mxa_₂F₁(n+1,n+α+1,2n+α+β+2,x))
stieltjesjacobimoment(α::Real,β::Real,n::Int,z) =
    (x = 2/(1-z);HypergeometricFunctions.mxa_₂F₁(n+1,n+α+1,2n+α+β+2,x))/2
stieltjesjacobimoment(α::Real,β::Real,z) = stieltjesjacobimoment(α,β,0,z)

function matanh(z)
    1/2*(log(1+z) - log(-1+z))
end

function matanh_p(z) # limit from above for z > 1
    1/2*(log(1+0im+z) - log(-1+0im+z)) + 1im*pi
end

function matanh_m(z) # limit from below for z > 1
    1/2*(log(1+0im+z) - log(-1+0im+z))
end

function legendrestieltjes(z)
    if real(z) <= 0
        return 1im/(4*pi)*(log(-complex(1)-z)-log(complex(1)-z))
    else
        -legendrestieltjes(-z)
    end
end

function legendrestieltjes_pos(z)
    return 1im/(4*pi)*(log(1+z) - 1im*pi - log(1-z))
end

function legendrestieltjes_neg(z)
    return 1im/(4*pi)*(log(1+z) + 1im*pi - log(1-z))
end

function MPSeedUnmapped(d,z)
    dp = (1 + sqrt(d))^2
    dm = (1 - sqrt(d))^2
    (1 - d - z + sqrt(z - dp)*sqrt(z - dm))/(4im*d*pi*z)
end

function MPSeedUnmappedPos(d,z)
    dp = (1 + sqrt(d))^2
    dm = (1 - sqrt(d))^2
    ((d -1 - z) + 1im*sqrt(dp - z)*sqrt(z - dm))/(4im*pi*z)
end

function MPSeedUnmappedNeg(d,z)
    dp = (1 + sqrt(d))^2
    dm = (1 - sqrt(d))^2
    ((d -1 - z) - 1im*sqrt(dp - z)*sqrt(z - dm))/(4im*pi*z)
end

function MPSeed(d)
    function f(z)
        Z = 2*sqrt(d)*z + 1 + d
        MPSeedUnmapped(d,Z)*2*sqrt(d)
    end
    f
end

function MPSeedPos(d)
    function f(z)
        Z = 2*sqrt(d)*z + 1 + d
        MPSeedUnmappedPos(d,Z)*2*sqrt(d)
    end
    f
end

function MPSeedNeg(d)
    function f(z)
        Z = 2*sqrt(d)*z + 1 + d
        MPSeedUnmappedNeg(d,Z)*2*sqrt(d)
    end
    f
end

function JacobiSeed(α,β)
    if α == 0.0 && β == 0.0
        return z -> legendrestieltjes(z)
    elseif α == 0.5 && β == 0.0 # log singularity at z = -1
        out = JacobiSeed(β,α)
        return z -> -out(-z)
    elseif α == 0.0 && β == 0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ?  3im/(8π)*(-2 + log(4)) + 3/(4*sqrt(2))*2*sqrt(2)*legendrestieltjes(z) :
         1/(8im*pi)*(6 - 3*sqrt(2)*sqrt(1+z + 0im)*matanh(sqrt(1+z+ 0im)/sqrt(2)))
    elseif α == -0.5 && β == 0.0
        out = JacobiSeed(β,α)
        z -> -out(-z)
    elseif α == 0.0 && β == -0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ? 1im*log(2.0)/(4π) + 0.5*legendrestieltjes(z)  :
        -1/(4im*sqrt(2)*pi)*1/sqrt(1+z +0im)*log((z-1)/(3 + complex(z) - 2*sqrt(2)*sqrt(1 + z |> complex)))
    else
        return z -> 1im/(2*pi)*stieltjesjacobimoment(α,β,z)
    end
end

function JacobiSeedPos(α,β)
    if α == 0.0 && β == 0.0
        return z -> legendrestieltjes_pos(z)
    elseif α == 0.5 && β == 0.0 # log singularity at z = -1
        out = JacobiSeedNeg(β,α)
        return z -> -out(-z)
    elseif α == 0.0 && β == 0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ?  3im/(8π)*(-2 + log(4)) + 3/(4*sqrt(2))*2*sqrt(2)*legendrestieltjes_pos(z) :
         1/(8im*pi)*(6 - 3*sqrt(2)*sqrt(1+z + 0im)*matanh_m(sqrt(1+z+ 0im)/sqrt(2)))
    elseif α == -0.5 && β == 0.0
        out = JacobiSeedNeg(β,α)
        return z -> -out(-z)
    elseif α == 0.0 && β == -0.5
        return z -> abs(z - 1) < 1e-14 ? 1im*log(2.0)/(4π) + 0.5*legendrestieltjes_pos(z)  :
        -1/(4im*sqrt(2)*pi)*1/sqrt(1+z)*(log(-(z-1)/(3 + z - 2*sqrt(2)*sqrt(1 + z))) - 1im*pi )
    else
        return z -> 1im/(2*pi)*stieltjesjacobimoment(α,β,z + 1im*eps())
    end
end

function JacobiSeedNeg(α,β)
    if α == 0.0 && β == 0.0
        return z -> legendrestieltjes_neg(z)
    elseif α == 0.5 && β == 0.0 # log singularity at z = -1
        out = JacobiSeedPos(β,α)
        z -> -out(-z)
    elseif α == 0.0 && β == 0.5 # log singularity at z = 1
        return z -> abs(z - 1) < 1e-14 ?  3im/(8π)*(-2 + log(4)) + 3/(4*sqrt(2))*2*sqrt(2)*legendrestieltjes_neg(z) :
        1/(8im*pi)*(6 - 3*sqrt(2)*sqrt(1+z)*(matanh_m(sqrt(1+z)/sqrt(2)) + 1im*pi))
    elseif α == -0.5 && β == 0.0
        out = JacobiSeedPos(β,α)
        z -> -out(-z)
    elseif α == 0.0 && β == -0.5
        return z -> abs(z - 1) < 1e-14 ? 1im*log(2.0)/(4π) + 0.5*legendrestieltjes_neg(z)  :
        -1/(4im*sqrt(2)*pi)*1/sqrt(1+z)*(log(-(z-1)/(3 + z - 2*sqrt(2)*sqrt(1 + z))) + 1im*pi )
    else
        return z -> 1im/(2*pi)*stieltjesjacobimoment(α,β,z - 1im*eps())
    end
end

function dist(z,n) # check if inside Bernstein ellipse that tends to
    # [-1,1] as n -> ∞
    ρ = 1 + 0.5/n
    a = (ρ+1/ρ)/2.0
    b = (ρ-1/ρ)/2.0
    if real(z)^2/a^2 + imag(z)^2/b^2 <= 1
        return 1
    else
        return 0
    end
end

@memoize function cauchy(a,b,seed,n,z::Number)
    if  dist(z,n) == 0   # the criterion for changing.
        m = 16;
        err = 1
        while err > 1e-15
            m *= 2
            c = fill(0.0im,m)
            c[1] = 1.0/(2im*pi)
            ldiv!(jacobi(a,b,m-1) - complex(z)*I,c)
            err = maximum(abs.(c[end-3:end]))
        end
        if m < n + 1
            append!(c,zeros(n + 10 - m))
        end
    else
        c = fill(0.0im,n+3)
        c[1] = seed(z);
        Z = complex(z)
        c[2] = Z*c[1] - a(0)*c[1] + 1/(2im*pi)
        c[2] = c[2]/b(0)
        for j = 1:n-1 # compute c_n
            c[j+2] = Z*c[j+1] - a(j)*c[j+1] - b(j-1)*c[j]
            c[j+2] /= b(j)
        end
    end
    c[1:n+1] 
end

function cauchy(a,b,seed,n,z::Vector)  # vectorize!
    A = zeros(ComplexF64,length(z),n+1)
    for i = 1:length(z)
        A[i,:] = cauchy(a,b,seed,n,z[i])
    end
    A
    #vcat(map(zz -> cauchy(a,b,seed,n,zz) |> transpose, z)...)
end

#### General recurrence coefficients ####
#### Note that this process is       ####
#### inherently unstable for         ####
#### unbounded intervals             ####
mutable struct RecCoef
    const ΛW::Function
    const NN::Function
    J::SymTridiagonal
    function RecCoef(ΛW,NN)
        new(ΛW,NN,SymTridiagonal([1.0,1.0],[0.0]))
    end
end

function _update!(RC::RecCoef,n::Int64)
    N = RC.NN(n+1)
    Λ, W = RC.ΛW(N)
    RC.J = lancz(n+1,Λ,W)
end

function _a(RC::RecCoef,n)
    if size(RC.J,1) < n + 1
        _update!(RC,2n)
    end
    RC.J[n+1,n+1]
end

function _b(RC::RecCoef,n)
    if size(RC.J,1) < n + 2
        _update!(RC,2n+1)
    end
    RC.J[n+1,n+2]
end

