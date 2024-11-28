# % LANCZOS Lanczos algorithm. --- https://www.cs.purdue.edu/archives/2002/wxg/codes/lanczos.m
# %
# %    Given the discrete inner product whose nodes are contained 
# %    in the first column, and whose weights are contained in the 
# %    second column, of the nx2 array xw, the call ab=LANCZOS(n,xw)
# %    generates the first n recurrence coefficients ab of the 
# %    corresponding discrete orthogonal polynomials. The n alpha-
# %    coefficients are stored in the first column, the n beta-
# %    coefficients in the second column, of the nx2 array ab.
# %
# %    The script is adapted from the routine RKPW in
# %    W.B. Gragg and W.J. Harrod, ``The numerically stable 
# %    reconstruction of Jacobi matrices from spectral data'', 
# %    Numer. Math. 44 (1984), 317-335.
# %
function lancz(N,x,w)
    Ncap = size(x,1);
    if N <= 0 || N > Ncap
        @error "lancz: N out of range"
    end
    p0 = copy(x)
    p1 = zeros(Ncap)
    p1[1] = w[1];
    for n = 1:Ncap-1
        pn = w[n+1]
        gam = 1
        sig=0
        t=0
        xlam = x[n+1];
        for k = 1:n+1
            rho = p1[k] + pn
            tmp = gam*rho
            tsig = sig
            if rho <= 0
                gam = 1
                sig = 0
            else
                gam = p1[k]/rho;
                sig = pn/rho;
            end
            tk = sig*(p0[k]-xlam) - gam*t
            p0[k] -= (tk-t)
            t = tk
            if sig <= 0
                pn = tsig*p1[k]
            else
            pn=(t^2)/sig;
            end
            tsig = sig
            p1[k] = tmp;
        end
    end
    SymTridiagonal(p0[1:N], sqrt.(p1[2:N]))
end

dot(f::BasisExpansion,g::BasisExpansion) = Base.sum(Multiplication(f)*Base.conj(g))

function dot(f::BasisExpansion{T},g::BasisExpansion{T}) where T <: DirectSum
    inner_prod = 0
    for i=1:length(f)
        inner_prod += dot(f[i],g[i])
    end
    return inner_prod
end

function sumdot(f::BasisExpansion,g::BasisExpansion)
    return dot(f,g)
end

#if f=[a,b], g=[c,d], this does <a,c> + <a,d> + <b,c> + <b,d>
function sumdot(f::BasisExpansion{T},g::BasisExpansion{T}) where T <: DirectSum
    inner_prod = 0
    for i=1:length(f)
        for j=1:length(g)
            inner_prod += dot(f[i],g[j])
        end
    end
    return inner_prod
end

function sumdot(f::BasisExpansion,g::BasisExpansion{T}) where T <: DirectSum
    inner_prod = 0
    for j=1:length(g)
        inner_prod += dot(f,g[j])
    end
    return inner_prod
end

function sumdot(f::BasisExpansion{T},g::BasisExpansion) where T <: DirectSum
    inner_prod = 0
    for i=1:length(f)
        inner_prod += dot(f[i],g)
    end
    return inner_prod
end

function sumdot(v1::Vector{T},v2::Vector{T}) where T <: BasisExpansion
    inner_prod = 0
    for j=1:length(v2)
        inner_prod += sumdot(v1[j],v2[j])
    end
    return inner_prod
end

function sumdot(v1::Vector,v2::Vector)
    inner_prod = 0
    for j=1:length(v2)
        inner_prod += sumdot(v1[j],v2[j])
    end
    return inner_prod
end