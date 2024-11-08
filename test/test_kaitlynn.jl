using OperatorApproximation

##################################################################################################################
#Nonoscillatory tests

#Function interpolation
N = 200;
α = 0;
gd = RationalRealAxis()
RationalRealAxis <: GridAxis
OscRational <: Basis
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp,N)
x_test = 0.1
abs(f(x_test) - ff(x_test)) < 1e-10
abs(f(x_test) - ff(x_test))

#Multplication
M = Multiplication(ff)
g = M*M*ff
abs(g(x_test)-ff(x_test)^3) < 1e-10
abs(g(x_test)-ff(x_test)^3)

#Derivative
D = Derivative()
h = D*ff
df = x-> -2*x*exp(-x^2)
abs(df(x_test) - h(x_test)) < 1e-10
abs(df(x_test) - h(x_test))

#Alternative way to test derivative
Op  = D*ff.basis
dc = Matrix(Op,N+2,N)*ff.c
test = BasisExpansion(ff.basis, dc)
abs(test(x_test)-df(x_test)) < 1e-10
abs(test(x_test)-df(x_test))


########################################################################################################################
#Oscillatory tests

#Function interpolation
N = 200;
α = 2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp,N)
x_test = 0.1
abs(f(x_test)*exp(1im*α*x_test) - ff(x_test)) < 1e-10
abs(f(x_test)*exp(1im*α*x_test) - ff(x_test))

#Multplication
M = Multiplication(ff)
g = M*M*ff
abs(g(x_test)*exp(1im*α*x_test)^2-ff(x_test)^3) < 1e-10
abs(g(x_test)*exp(1im*α*x_test)^2-ff(x_test)^3)

#Derivative
D = Derivative()
h = D*ff
df = x-> -2*x*exp(-x^2)*exp(1im*α*x) + exp(-x^2)*1im*α*exp(1im*α*x)
abs(df(x_test) - h(x_test)) < 1e-10
abs(df(x_test) - h(x_test))

##########################################################################################################################
#Cauchy tests

###########################################
#C^+ tests

#Computes the "true" C+ integral using quadrature so that we don't have to keep comparing with Mathematica
function CPquad(f,a,z)
    R = 20
    ff = x -> f(x)*exp(1im*x*a)/(2im*pi)
    F = x -> f(x)*exp(1im*x*a)/(2im*pi)*1/(x-z)
    endpts = [-R, z-1, z-1im, z+1, R]
    s = 0
    for i = 1:4
        gd = JacobiMappedInterval(endpts[i],endpts[i+1],0,0)
        sp = Jacobi(0,0,gd)
        ff = BasisExpansion(F,sp)
        s = s + sum(ff)
    end
    s
end

#Nonoscillatory tests
α = 0;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
C_plus = CauchyOperator(1)
x_test = 0.145;

##Even number of coefficients
N = 200;
ff = BasisExpansion(f,sp,N)
(C_plus*ff)(x_test) - CPquad(f,α,x_test)

##Odd number of coefficients
N = 201;
ff = BasisExpansion(f,sp,N)
(C_plus*ff)(x_test) - CPquad(f,α,x_test)

#Oscillatory test (α > 0)
α = 2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
C_plus = CauchyOperator(1)
x_test = 0.145;

##Even number of coefficients
N = 200;
ff = BasisExpansion(f,sp,N)
(C_plus*ff)(x_test) - CPquad(f,α,x_test)

##Odd number of coefficients
N = 201;
ff = BasisExpansion(f,sp,N)
(C_plus*ff)(x_test) - CPquad(f,α,x_test)

#Oscillatory test (α < 0)
α = -2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
C_plus = CauchyOperator(1)
x_test = 0.145;

##Even number of coefficients
N = 200;
ff = BasisExpansion(f,sp,N)
(C_plus*ff)(x_test) - CPquad(f,α,x_test)

##Odd number of coefficients
N = 201;
ff = BasisExpansion(f,sp,N)
(C_plus*ff)(x_test) - CPquad(f,α,x_test)

###########################################
#C^- tests

#Computes the "true" C+ integral using quadrature so that we don't have to keep comparing with Mathematica
function CMquad(f,a,z)
    R = 20
    ff = x -> f(x)*exp(1im*x*a)/(2im*pi)
    F = x -> f(x)*exp(1im*x*a)/(2im*pi)*1/(x-z)
    endpts = [-R, z-1, z+1im, z+1, R]
    s = 0
    for i = 1:4
        gd = JacobiMappedInterval(endpts[i],endpts[i+1],0,0)
        sp = Jacobi(0,0,gd)
        ff = BasisExpansion(F,sp)
        s = s + sum(ff)
    end
    s
end

#Nonoscillatory tests
α = 0;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
C_minus = CauchyOperator(-1)
x_test = 0.145;

##Even number of coefficients
N = 200;
ff = BasisExpansion(f,sp,N)
(C_minus*ff)(x_test) - CPquad(f,α,x_test) #need to use CPquad to test here b/c that is what is used for α>=0

##Odd number of coefficients
N = 201;
ff = BasisExpansion(f,sp,N)
(C_minus*ff)(x_test) - CPquad(f,α,x_test) #need to use CPquad to test here b/c that is what is used for α>=0

#Oscillatory test (α > 0)
α = 2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
C_minus = CauchyOperator(-1)
x_test = 0.145;

##Even number of coefficients
N = 200;
ff = BasisExpansion(f,sp,N)
(C_minus*ff)(x_test) - CMquad(f,α,x_test)

##Odd number of coefficients
N = 201;
ff = BasisExpansion(f,sp,N)
(C_minus*ff)(x_test) - CMquad(f,α,x_test)

#Oscillatory test (α < 0)
α = -2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
C_minus = CauchyOperator(-1)
x_test = 0.145;

##Even number of coefficients
N = 200;
ff = BasisExpansion(f,sp,N)
(C_minus*ff)(x_test) - CMquad(f,α,x_test)

##Odd number of coefficients
N = 201;
ff = BasisExpansion(f,sp,N)
(C_minus*ff)(x_test) - CMquad(f,α,x_test)

##############################################################################
#Integral of an OscRational basis expansion tests
α = -2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2);
ff = BasisExpansion(f,sp)
int_test = Base.sum(ff)
int_true = sqrt(π)*exp(-(abs(α)^2)/4)
abs(int_test-int_true)

f2 = x -> exp(-2x^2);
ff2 = BasisExpansion(f2,sp)
int_test2 = Base.sum(ff2)
int_true2 = sqrt(π/2)*exp(-(abs(α)^2)/8)
abs(int_test2 - int_true2)

###############################################################################
#Multiplication on bases with two different alpha values

α = 0;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp)
M = Multiplication(ff)

α2 = 1;
sp2 = OscRational(gd,α2)
gg = BasisExpansion(f,sp2)

x_test = 0.356
h = M*gg
abs(h(x_test)-ff(x_test)*gg(x_test)) < 1e-10
abs(h(x_test)-ff(x_test)*gg(x_test))

#################################################################################
#Conjugate of a basis expansion
α = 0;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> x/(x+1im)
ff = BasisExpansion(f,sp)
x_test = 5.67
abs(f(x_test)-ff(x_test))

ff2 = Base.conj(ff)
conj_f = x -> x/(x-1im)
abs(ff2(x_test) - conj_f(x_test))

###################################################################################
#Test dot product of BasisExpansion{OscRational}
α = 2;
gd = RationalRealAxis()
sp = OscRational(gd,α)
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp)

α = -3;
sp = OscRational(gd,α)
g = x -> exp(-2x^2)
gg = BasisExpansion(g,sp)

dot(ff,gg)
sqrt(pi/3)*exp(-25/12) #have to take into account exponential factors in integrand as well
abs(dot(ff,gg) - sqrt(pi/3)*exp(-25/12))

#####################################################################################
#Test dot product of BasisExpansion{DirectSum}
α = 2;
gd = RationalRealAxis()
sp = OscRational(gd,α)
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp)

α = -3;
sp = OscRational(gd,α)
g = x -> exp(-2x^2)
gg = BasisExpansion(g,sp)

h1 = ff ⊕ gg
h2 = gg ⊕ ff

dot(h1,h2)
2*sqrt(pi/3)*exp(-25/12) #\int_{-∞}^{+∞} ff*conj(gg) dx = \int_{-∞}^{+∞} gg*conj(ff) dx (hence the times 2 from answer above)
abs(dot(h1,h2) - 2*sqrt(pi/3)*exp(-25/12))