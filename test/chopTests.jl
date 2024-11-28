using OperatorApproximation

import Base: copy

function copy(f::BasisExpansion)
    BasisExpansion(f.basis,copy(f.c))
end

#comptue ψ (sign of L and σ determines if 1,2,+,-)
function ψ(V,L,k,σ)
    if L > 0
        A = 0.0;
        B = L;
    else
        A = L;
        B = 0.0;
    end
    
    grid = ChebyshevMappedInterval(A,B)
    sp = Ultraspherical(0.0,grid)
    D = Derivative()
    gv = GridValues(grid)
    E = Conversion(gv)
    M = Multiplication(V)
    Op = -E*D^2 - (2*1im*σ*k)*E*D + M*E
    bdryE = FixedGridValues([L],grid) |> Conversion
    L1 = (bdryE ⊘ (bdryE*D) ⊘ Op)*sp
    N_sol = (L1 \ [[0];[0]; x -> -V(x)])[1]
    dN_sol = D*N_sol
    
    ψn = x-> (N_sol(x)+1)*exp(1im*k*σ*x)
    dψn = x-> (N_sol(x)+1)*(1im*σ*k)*exp(1im*k*σ*x) + exp(1im*k*σ*x)*dN_sol(x)
    
    return ψn,dψn
end

#compute ϕ (solutions to homogeneous problem of ψ bounded at ±∞ depending on sign of L)
function Φ(V,L,k,f)
    if L > 0
        A = 0.0;
        B = L;
    else
        A = L;
        B = 0.0;
    end
    
    grid = ChebyshevMappedInterval(A,B)
    sp = Ultraspherical(0.0,grid)
    D = Derivative()
    gv = GridValues(grid)
    E = Conversion(gv)
    M = Multiplication(x-> V(x)-k^2)
    Op = -E*D^2 + M*E
    bdryE = FixedGridValues([L],grid) |> Conversion
    L1 = (bdryE ⊘ (bdryE*D) ⊘ Op)*sp
    Φn = (L1 \ [[0];[0];f])[1]
    dΦn = D*Φn
    
    return Φn, dΦn
end

#computes C(k) (forward transform without constant out front)
function C(V,L,k,f)
    x = 0.0
    ψ1_minus, dψ1_minus = ψ(V,-L,k,-1.0)
    ψ2_plus, dψ2_plus = ψ(V,L,k,1.0)
    Φ_plus, dΦ_plus = Φ(V,L,k,f) 
    Φ_minus, dΦ_minus = Φ(V,-L,k,f) 
    out = -1*((Φ_plus(x) - Φ_minus(x))*dψ1_minus(x) - ψ1_minus(x)*(dΦ_plus(x) - dΦ_minus(x))) 
end

#compute the Wronskian
function W(a,b)
    z = 0.0
    a[1](z)*b[2](z) - a[2](z)*b[1](z)
end

#computes the terms of the scattering matrix
function S(V,L,k) 
    ψ1p =  ψ(V,-L,k,1.0)
    ψ1m =  ψ(V,-L,k,-1.0)
    ψ2p =  ψ(V,L,k,1.0)
    ψ2m =  ψ(V,L,k,-1.0)
    a11 = W(ψ1m,ψ2p)/W(ψ2m,ψ2p)
    a12 = W(ψ1m,ψ2m)/W(ψ2p,ψ2m)
    a21 = W(ψ1p,ψ2p)/W(ψ2m,ψ2p)
    a22 = W(ψ1p,ψ2m)/W(ψ2p,ψ2m)
    [a11 a12; a21 a22]
end

#computes reflection coefficient, ρ1, for x<0
function r(k)
    if abs(k) ≈ 0.0
        return -1.0+0.0im
    elseif abs(k) > 100
        return 0.0*1im
    else
        Scat_Mat = S(V,L,k) 
        return -Scat_Mat[2,1]/Scat_Mat[1,1]
    end 
end

#computes reflection coefficient, ρ2, for x<0
function rbar(k)
    if abs(k) ≈ 0.0
        return -1.0+0.0im
    elseif abs(k) > 100
        return 0.0*1im
    else
        Scat_Mat = S(V,L,k)
        return -Scat_Mat[1,2]/Scat_Mat[2,2]
    end 
end

#zai = zero at infinity [had to track this down from ApproxFunRational]
function zai(f)
    return x -> isnan(f(x)) ? 0.0im : complex(f(x))
end

###################################################################################################################################
n = 200
V = x -> (sech(x))^2;
L = 20;
f = x -> exp(-x^2);

gd = RationalRealAxis()
sp = OscRational(gd,0.0) #using OscRational here instead of OscLaurent from ApproxFun
rfun = BasisExpansion(k -> r(k),sp,n) |> Base.chop #using BasisExpansion instead of Fun 
rbarfun = BasisExpansion(k -> rbar(k),sp,n) |> Base.chop
# Cosc = BasisExpansion(zai(k -> C(V,L,k,f)),sp,800) |> Base.chop #kind of need the extra coefficients on this one to effectively capture it, but slows GMRES down
Cosc = BasisExpansion(zai(k -> C(V,L,k,f)),sp,n) |> Base.chop

#############################################################################################################
x = -1

α = -2*x;
sp = OscRational(rfun.basis.GD,α)
ρ1 = BasisExpansion(sp,rfun.c)

sp = OscRational(rbarfun.basis.GD,-α)
ρ2 = BasisExpansion(sp,rbarfun.c)

M = Multiplication(-1*ρ1)
h = [M*Cosc ⊕ 0*Cosc,Cosc ⊕ 0*Cosc]
# h = [M*Cosc, Cosc]

###################################################################################################################
maximum(abs.(h[1][2].c)) < 1e-15

hSimp = simp(h)
hSimp[1][1].c
maximum(abs.(hSimp[1][2].c))
length(hSimp[1])

test1 = Base.chop(h)
test2 = combine(test1)
test3 = Base.chop(test2)
####################################################################################################################

#computes f(x) for x<0
function computef(x,tol,rfun,rbarfun,Cosc)
    α = -2*x;
    sp = OscRational(rfun.basis.GD,α)
    ρ1 = BasisExpansion(sp,rfun.c)

    sp = OscRational(rbarfun.basis.GD,-α)
    ρ2 = BasisExpansion(sp,rbarfun.c)

    M = Multiplication(-1*ρ1)
    h = [M*Cosc ⊕ 0*Cosc,Cosc ⊕ 0*Cosc]
    # h = Any[M*Cosc,Cosc]

    h = simp(h)
       
    𝓒⁺ = CauchyOperator(1)
    𝓒⁻ = CauchyOperator(-1)

    M2 = Multiplication(ρ2)
    function u(x)
        copied = copy.(x)
        out = Vector{BasisExpansion}(undef,2)
        out[1] = copied[1]
        out[2] = copied[2]
        # out = copy.(x)
        # out = x
        # out[1] = out[1] ⊕ M*(𝓒⁺*x[2])
        # out[2] = out[2] ⊕ M2*(𝓒⁻*x[1]) 
        out[1] = out[1] + M*(𝓒⁺*x[2])
        out[2] = out[2] + M2*(𝓒⁻*x[1]) 
        return out
    end

    # out = GMRES(u,h,⋅,tol,30,x -> simp(x)) #should I keep the last entry as x-> simp(x)
    out = GMRES(u,h,sumdot,tol,30,x -> simp(x)) #something has to be done about this simp() function...
   #  sol = +([out[2][i]*out[1][i] for i=1:length(out[2])]...)
        
   # out = sum(sol)/(2*π)
   # -1*out[1]+out[2]
end

function GMRES(A,b,inner,tol,n,cond)
    println("I am in GMRES!")
    nom = a -> sqrt(abs(inner(a,a)))
    H = zeros(Complex{Float64},n+1,n)
    bnorm = nom(b)

    x = 0.
    conv_history = []
    Q = [(1.0/bnorm)*b]

    for i = 1:n
       #tic()
       #println("Operator application: ")
       v = A(Q[i])
       #toc()
       #tic()
       #println("Inner products: ")
       for j = 1:i
           H[j,i] = inner(Q[j],v)
           v = cond(v - H[j,i]*Q[j])
       end
       v = cond(v)
       #println("Assembling Q:")
       H[i+1,i] = nom(v)
       Q = vcat(Q,[copy((1.0/H[i+1,i])*v)])
       #print("Arnoldi: ")
       #toc()
       #return v
       if i > 0
           # Solve H[1:i+1,1:i]*x = bnorm*e_1, using least squares
           # TODO: Implement Givens rotations
           rhs = zeros(Float64,i+1)
           rhs[1] = bnorm
           x = H[1:i+1,1:i]\rhs
           res = norm(H[1:i+1,1:i]*x-rhs)
           conv_history = vcat(conv_history,[i,res])
           print("iteration = ")
           print(i)
           print(", residual = ")
           println(res)
           if res < tol
               return  [Q,x,conv_history]
           end
       end
    end
    println("GMRES did not terminate")
    return [Q,x,conv_history]
end
###################################################################################################
out = computef(-4,1e-5,rfun,rbarfun,Cosc)

###################################################################################################
#scalar test RHP
# phiPlus = phiMinus*(1+exp(-x^2))

using OperatorApproximation

import Base: copy

function copy(f::BasisExpansion)
    BasisExpansion(f.basis,copy(f.c))
end

gd = RationalRealAxis()
sp = OscRational(gd,0)
G = x -> 1+0.01*exp(-x^2)
# G = x -> 1+(1/x^2)
𝓒⁺ = CauchyOperator(1)
𝓒⁻ = CauchyOperator(-1)
GMinus1 = BasisExpansion(x->G(x)-1,sp)
GMinus1 = chop(GMinus1)
M = Multiplication(GMinus1)
Sop = u-> u-M*𝓒⁻*u

sol = GMRES(Sop,GMinus1,sumdot,1e-10,30,chop)

u = sum(sol[2].*sol[1][1:end-1])
(𝓒⁺*u)(0)+1
((𝓒⁻*u)(0)+1)*G(0)

G1 = BasisExpansion(x->(1/G(x))-1,sp)
G1 = chop(G1)
M1 = Multiplication(G1)
G1.c
test = chop(G1)
test.c
S1 = u-> u-M1*𝓒⁻*u
S2 = u-> S1(Sop(u))
S2(GMinus1)
sol1 = GMRES(S2,S1(GMinus1),sumdot,1e-10,30,chop)

#########################################################################

# gd = RationalRealAxis();
# sp = OscRational(gd,1.0);
# w = BasisExpansion(x -> 2im/(x + 1im)^2, sp);
# f = BasisExpansion(x -> exp(-x^2),sp);
# M1 = Multiplication(w)
# M = Multiplication(w)*sp;
# (M*f)(1.0) - f(1.0)*w(1.0)
# test = M*f
# test(1.0)
# w(1.0)*f(1.0)