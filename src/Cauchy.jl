function dist(z,n) # check if inside Bernstein ellipse that tends to
    # [-1,1] as n -> ∞
    ρ = 1 + 1/n
    a = (ρ+1/ρ)/2.0
    b = (ρ-1/ρ)/2.0
    if real(z)^2/a^2 + imag(z)^2/b^2 <= 1
        return 1
    else
        return 0
    end
end


# This is all code to do adaptive QR with a tridiagonal matrix

mutable struct QuadDiagonal  # length(diag) >= 3
   sub::Vector
   diag::Vector
   sup::Vector
   supsup::Vector
end

function Array(A::QuadDiagonal)
    b = (-1 => A.sub, 0 => A.diag, 1 => A.sup, 2 => A.supsup)
    diagm(b...)
end

function toBandedMatrix(A::QuadDiagonal)
    b = (-1 => A.sub, 0 => A.diag, 1 => A.sup, 2 => A.supsup)
    BandedMatrix(b,(length(A.diag),length(A.diag)))
end

function extend!(A::QuadDiagonal,v::Vector)
    append!(A.sub,v[1])
    append!(A.diag,v[2])
    append!(A.sup,v[3])
    append!(A.supsup,v[4])
end

function complex(A::QuadDiagonal)
    QuadDiagonal(A.sub |> complex, A.diag |> complex, A.sup |> complex, A.supsup |> complex) 
end

function complex!(A::QuadDiagonal)
    A.sub = A.sub |> complex
    A.diag = A.diag |> complex
    A.sup = A.sup |> complex
    A.supsup = A.supsup |> complex
end

function givens!(A::QuadDiagonal,j)
    if j > length(A.diag) -1
        @error "Matrix too small"
        return [1,0]
    end
    a = A.diag[j]
    b = A.sub[j]
    nom = sqrt(abs2(a)+abs2(b))
    c = a/nom; s = b/nom
    A.diag[j] = nom
    A.sub[j] = 0.0
    if j < n-1
        if A.supsup[j] != 0.0
        @warn "Structure issue, factorization not correct"
        end
        R = [conj(c) conj(s); -s c]*[A.sup[j] 0.0; A.diag[j+1] A.sup[j+1]]
        A.sup[j] = R[1,1]
        A.supsup[j] = R[1,2]
        A.diag[j+1] = R[2,1]
        A.sup[j+1] = R[2,2]
    else
        R = [conj(c) conj(s); -s c]*[A.sup[j] ; A.diag[j+1]]
        A.sup[j] = R[1]
        A.diag[j+1] = R[2]
    end
    [c,s]
end

function givens!(A::QuadDiagonal,v::Vector,j)
    n = length(A.diag)
    if j > n -1 || length(v) != n
        @error "Matrix wrong size"
        return A,v
    end
    a = A.diag[j]
    b = A.sub[j]
    nom = sqrt(abs2(a)+abs2(b))
    c = a/nom; s = b/nom
    A.diag[j] = nom
    A.sub[j] = 0.0
    if j < n-1
        if A.supsup[j] != 0.0
            @warn "Structure issue, factorization not correct"
        end
        R = [conj(c) conj(s); -s c]*[A.sup[j] 0.0 v[j]; A.diag[j+1] A.sup[j+1] v[j+1]]
        A.sup[j] = R[1,1]
        A.supsup[j] = R[1,2]
        A.diag[j+1] = R[2,1]
        A.sup[j+1] = R[2,2]
        v[j] = R[1,3]
        v[j+1] = R[2,3]
    else
        R = [conj(c) conj(s); -s c]*[A.sup[j] v[j]; A.diag[j+1] v[j+1]]
        A.sup[j] = R[1,1]
        A.diag[j+1] = R[2,1]
        v[j] = R[1,2]
        v[j+1] = R[2,2]
    end
    [c,s]
end

function QuadDiagonal(J::SymTridiagonal)
   QuadDiagonal(diag(J,-1) |> complex, diag(J) |> complex, diag(J,1) |> complex, fill(0.0im,length(diag(J))-2)) 
end

function cauchy_off(a,b,nn,z)
    n = 3
    A = Jacobi(a,b,n) - z*I |> QuadDiagonal;
    v = fill(0.0im,n+1)
    v[1] = 1.0/(2im*pi)
    cs = [0im,0im]
    for i = 1:n
        cs = givens!(A,v,i)
    end
    n += 1
    while abs(v[end]) > 1e-16 && n < 10000
        c = cs[1]; s = cs[2]
        #R = [conj(c) conj(s); -s c]*[0im, b(n)]
        #R = b(n)*[conj(s), c]
        extend!(A,[b(n),a(n+1)-z, b(n)*c , b(n)*conj(s)])
        append!(v,0.0)
        cs = givens!(A,v,n)
        n += 1
    end
    pad(toBandedMatrix(A)\v,nn+1)
end

function cauchy(a,b,seed,n,z::Number)
    if  dist(z,n) == 0   # the criterion for changing.
        # Non-adaptive method
        #m = 2; #over sampling, should be done adaptively
        #v = fill(0.0im,m*n+1)
        #v[1] = 1.0/(2im*pi)
        #c = ((Jacobi(a,b,m*n) - z*I)\v)
        # Adaptive method
        #println("eval off")
        c = cauchy_off(a,b,n,z)
    else
        c = fill(0.0im,n+3)
        c[1] = seed(z);
        c[2] = z*c[1] - a(0)*c[1] + 1/(2im*pi)
        c[2] = c[2]/b(0)
        for j = 1:n-1 # compute c_n
            c[j+2] = z*c[j+1] - a(j)*c[j+1] - b(j-1)*c[j]
            c[j+2] /= b(j)
        end
    end
    c[1:n+1]
end

function cauchy(a,b,seed,n,z::Vector)
    vcat(map(zz -> cauchy(a,b,seed,n,zz) |> transpose, z)...)
end