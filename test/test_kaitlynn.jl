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

#Nonoscillatory test (Working properly)
N = 200;
α = 0;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp,N)
C_plus = CauchyOperator(1)
test1 = C_plus*ff
test1(0.145)

#Oscillatory test (α > 0)
N = 200;
α = 2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp,N)
C_plus = CauchyOperator(1)
test2 = C_plus*ff
test2(0.145)
test2(0.145) + ff(0.145)
test2(2)

#Oscillatory test (α < 0)
N = 201;
α = -2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp,N)
C_plus = CauchyOperator(1)
test3 = C_plus*ff
test3(0.145)
test3(2)