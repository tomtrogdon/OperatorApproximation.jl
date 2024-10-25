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