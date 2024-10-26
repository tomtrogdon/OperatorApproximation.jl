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

N = 200;
α = 0;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp,N)
Cp = CauchyOperator(1)
Cm = CauchyOperator(-1)
test1 = Cp*ff
@test abs(test1(0.145) - CPquad(f,α,0.145)) < 1e-10
test1 = Cm*ff
@test abs(test1(0.145) - CMquad(f,α,0.145)) < 1e-10

CPquad(f,α,0.145) - CMquad(f,α,0.145) - f(0.145)

α = 2;
gd = RationalRealAxis()
sp = OscRational(gd,α);
f = x -> exp(-x^2)
ff = BasisExpansion(f,sp,N)
Cp = CauchyOperator(1)
Cm = CauchyOperator(-1)
test1 = Cp*ff
@test abs(test1(0.145) - CPquad(f,α,0.145)) < 1e-10
test1 = Cm*ff
@test abs(test1(0.145) - CMquad(f,α,0.145)) < 1e-10