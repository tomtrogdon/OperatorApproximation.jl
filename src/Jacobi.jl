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

function Jacobi(a,b,n) # creates (n + 1) x (n+1) Jacobi matrix
   SymTridiagonal([a(i) for i in 0:n],[b(i) for i in 0:n-1])
end

function Interp_transform(a,b,n)
    E = Jacobi(a,b,n) |> eigen
    return E.values, E.vectors*(Diagonal(E.vectors[1,:]))
end

function Gauss_quad(a,b,n)
    E = Jacobi(a,b,n) |> eigen
    return E.values, abs2.(E.vectors[1,:])
end

function Interp_transform(f::Function,a,b,n)
    E = Jacobi(a,b,n) |> eigen
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