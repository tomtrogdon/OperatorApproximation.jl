using OperatorApproximation
gd = RationalRealAxis();
sp = OscRational(gd,0.0)
ff = BasisExpansion(x -> exp(-x^2), sp, 300)
spα = OscRational(gd,1.0)
ffα = BasisExpansion(x -> exp(-x^2), spα, 300)
fadd = ff + ffα;
Cop = CauchyOperator(1)*fadd.basis;
out = Cop*fadd

function CPquad(f,a,z)
    R = 10
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

CPquad(x -> exp(-x^2),1,3) + CPquad(x -> exp(-x^2),0,3) - out(3)




out(3)

out.basis[3]

  \\


Cop.L.Ops

OperatorApproximation.op_momtm(Cop)

OperatorApproximation.grabOps(ZeroOperator())


Cop[1,:][2]


OperatorApproximation.grabOps.(Cop)

typeof(Cop[1,1])

OperatorApproximation.grabOps(Cop[2,2])

Cop[2,2]

Cop.L.Ops[2,2].Ops



bb = DirectSum([fadd.basis, fadd.basis]);
bb.bases


CauchyOperator(1)*ffα.basis;

ffα.basis

Cop.range[2]

dC = CauchyOperator(1) ⊘ CauchyOperator(1)
dC*ffα.basis

dC.Ops

(AbstractZeroOperator()*ff.basis).range


b = DirectSum([OperatorApproximation.AnyBasis(), fadd.basis])

typeof(OperatorApproximation.AnyBasis()) <: OperatorApproximation.AnyBasis
typeof(fadd.basis) <: OperatorApproximation.AnyBasis

OperatorApproximation.AnyBasis() != ff.basis

b.bases

Cop2 = dC*fadd.basis
Cop2[2,2].range

Cop[2,2].L




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