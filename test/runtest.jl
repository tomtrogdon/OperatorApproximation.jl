using OperatorApproximation, SpecialFunctions, LinearAlgebra
using Test

@testset "OperatorApproximation.jl: Basic tests" begin
    R = 30;
    sp = Ultraspherical(0.0,ChebyshevMappedInterval(-R,R));
    f = BasisExpansion(sp,[0,1])
    M = Multiplication(f)
    g = M*M*f
    @test abs(g(.3) - f(.3)^3) < 1e-10

    R = 10;
    sp = Ultraspherical(0.0,ChebyshevMappedInterval(-R,R));
    f = BasisExpansion(sp,randn(10))
    sp = Ultraspherical(0.0,ChebyshevMappedInterval(-R,R));
    g = Conversion(sp)*f
    @test abs(g(1.0) - f(1.0)) < 1e-10

    sp = Ultraspherical(1.0,ChebyshevMappedInterval(-R,R));
    g = Conversion(sp)*f
    @test abs(g(1.0) - f(1.0)) < 1e-10

    sp = Ultraspherical(2.0,ChebyshevMappedInterval(-R,R));
    g = Conversion(sp)*f
    @test abs(g(1.0) - f(1.0)) < 1e-10

    sp = Ultraspherical(3.0,ChebyshevMappedInterval(-R,R));
    g = Conversion(sp)*f
    @test abs(g(1.0) - f(1.0)) < 1e-10

    f = x -> sin(x)
    ff = BasisExpansion(f,Ultraspherical(0.0,UltraInterval(0.0)),30)
    @test abs(ff(.3) - f(.3)) < 1e-10

    f = x -> sin(x)
    ff = BasisExpansion(f,Ultraspherical(1.0,UltraInterval(1.0)),30)
    @test abs(ff(.3) - f(.3)) < 1e-10

    gd = UltraMappedInterval(-10,10,0.0)
    f = x -> sin(x)
    ff = BasisExpansion(f,Ultraspherical(0.0,gd),100)
    @test abs(ff(3.0)-sin(3.0)) < 1e-10

    N = 10;
    sp = Fourier(PeriodicInterval());
    f = x -> sin(2*pi*x)
    ff = BasisExpansion(f,sp,N)
    @test abs(f(.1) - ff(.1)) < 1e-10

    N = 100
    f = x -> exp(-x^2)
    sp = Fourier(PeriodicMappedInterval(-10,10))
    ff = BasisExpansion(f,sp,N)
    @test abs(ff(.1) - f(.1)) < 1e-10

    f = x -> sin(x)
    ff = BasisExpansion(f,Laurent(PeriodicMappedCircle(1im,.5)))
    cff = CauchyTransform()*ff
    @test abs(cff(1.1im) - f(1.1im)) < 1e-10
end


@testset "OperatorApproximation.jl: Airy equation" begin
    # test collocation
    R = 30;
    sp = Ultraspherical(0.0,ChebyshevMappedInterval(-R,R));
    gv = GridValues(ChebyshevMappedInterval(-R,R));
    E = Conversion(gv);
    M = Multiplication(x -> x);
    D = Derivative();
    Op = E*D^2 - M*E
    lbdry = FixedGridValues([-R],ChebyshevMappedInterval(-R,R)) |> Conversion;
    rbdry = FixedGridValues([R],ChebyshevMappedInterval(-R,R)) |> Conversion;
    setbasis(sp)
    u = (lbdry ⊘ rbdry ⊘ Op)\[[airyai(-R)];[airyai(R)]; x->0]
    u = u[1]
    @test abs(u(0) - airyai(0)) < 1e-10

    # test sparse method
    R = 30;
    sp = Ultraspherical(0.0,UltraMappedInterval(-R,R,2.0)); #choose this grid domain
    # because it is inherited by the inferred spaces
    # so it needs to be compatible for a projection of the
    # right-hand side function
    sp2 = Ultraspherical(2.0,UltraMappedInterval(-R,R,2.0));
    M = Multiplication(x -> x);
    D = Derivative();
    Op = D^2 - Conversion(sp2)*M
    lbdry = FixedGridValues([-R],ChebyshevMappedInterval(-R,R)) |> Conversion;
    rbdry = FixedGridValues([R],ChebyshevMappedInterval(-R,R)) |> Conversion;
    setbasis(sp)
    u = (lbdry ⊘ rbdry ⊘ Op)\[[airyai(-R)];[airyai(R)]; x->0]
    u = u[1]
    @test abs(u(0) - airyai(0)) < 1e-10
end

@testset "OperatorApproximation.jl: Cauchy tests" begin
    # Note: This computes integrals
    #
    # \frac{1}{2 \pi i} \int_a^b \frac{f(s)}{s -z} w(s) d s
    #
    # where w(s) \geq 0 and
    #
    # \left| \int_a^b w(s) ds \right| = |b - a|
    #
    truth = 0.009366780921585525*im  # true integral wrt normalized weight function
    f = x -> 0.5*sin.(x)
    gd = JacobiMappedInterval(-1.0,1.0,0.5,0.5)
    sp = Jacobi(0.5,0.5,gd)
    ff = BasisExpansion(f,sp,100)
    @test abs(ff(.3) - f(.3)) < 1e-10
    cff = CauchyTransform()*ff
    z = 2.1
    @test abs(cff(z) - truth) < 1e-10

    # c = NIntegrate[Sqrt[4 - x^2], {x, -2, 2}]
    # truth = 1/(2 Pi I*c) 4*NIntegrate[Sin[x] Sqrt[4 - x^2]/(x - 2.1), {x, -2, 2}]
    truth = 0.21542074482463697*im  # true integral wrt normalized weight function
    f = x -> sin.(x)
    gd = JacobiMappedInterval(-2.0,2.0,0.5,0.5)
    sp = Jacobi(0.5,0.5,gd)
    ff = BasisExpansion(f,sp,100)
    @test abs(ff(.3) - f(.3)) < 1e-10
    cff = CauchyTransform()*ff
    z = 2.1
    display(cff(z)/truth)
    @test abs(cff(z)-truth) < 1e-10

    truth = 0.0808961206892101 -  0.1554484139468048*im
    f = x -> 0.5*sin.(x)*sqrt(2)
    gd = JacobiMappedInterval(0.0,1im + 1.0,0.5,0.5)
    sp = Jacobi(0.5,0.5,gd)
    ff = BasisExpansion(f,sp,100)
    @test abs(ff(.3im + .3) - f(.3im + .3)) < 1e-10
    cff = CauchyTransform()*ff
    @test abs(cff(.1) - truth) < 1e-10

    f = x -> 0.5*sin.(x)
    gd0 = JacobiMappedInterval(-1.0,1.0,0.0,0.0)
    sp = Jacobi(0.0,0.0,gd0)
    ff = BasisExpansion(f,sp,100)

    gd = DirectedLobattoMappedInterval(-1,1) 

    Bp = BoundaryValue(1,GridValues(gd))
    Cp = Bp*CauchyTransform()
    CCp = Cp*ff.basis

    Bm = BoundaryValue(-1,GridValues(gd))
    Cm = Bm*CauchyTransform()
    CCm = Cm*ff.basis

    E = Conversion(GridValues(gd))
    CE = E*ff.basis

    @test norm(Matrix(CCp,10,10) - Matrix(CCm,10,10) - Matrix(CE,10,10)) < 1e-10

end


@testset "OperatorApproximation.jl: Riemann-Hilbert test" begin
    gd = PeriodicCircle()
    sp = Laurent(gd)
    g = z -> 3 + z
    M = Multiplication(g)
    Cp = CauchyOperator(1)
    Cm = CauchyOperator(-1)
    setbasis(sp)
    u = \(Cp - M*Cm, z -> g(z) -1)
    m = CauchyTransform()*u
    @test abs(m(2im)) < 1e-10 && abs(m(.3) + 1 - g(.3)) < 1e-10

    g = x -> 1 + 0.5exp(-30x^2)
    f = x -> g(x) - 1
    gd1 = JacobiMappedInterval(-1,1.0,0.0,0.0)
    s1 = Jacobi(0.0,0.0,gd1)

    dgd1 = DirectedLobattoMappedInterval(-1,1)
    gv1 = GridValues(dgd1)

    Cp = BoundaryValue(1,gv1)*CauchyTransform() 
    Cm = BoundaryValue(-1,gv1)*CauchyTransform() 
    M = Multiplication(g)
    MCm = M*Cm
    S = Cp - MCm
    cS = S*s1
    u = \(cS,f,100)
    truth = 1.0220013279200346
    @test abs(truth - (CauchyTransform()*u)(1im) -1) < 1e-10

    g = x -> 1 + 0.5exp(-30x^2)
    f = x -> g(x) - 1
    gd1 = JacobiMappedInterval(-1,0,0.0,0.0)
    gd2 = JacobiMappedInterval(0,1,0.0,0.0)
    s1 = Jacobi(0.0,0.0,gd1)
    s2 = Jacobi(0.0,0.0,gd2)

    dgd1 = DirectedLobattoMappedInterval(-1,0)
    dgd2 = DirectedLobattoMappedInterval(0,1)
    gv1 = GridValues(dgd1)
    gv2 = GridValues(dgd2)

    r1 = BoundaryValue(1,gv1)*CauchyTransform() ⊞ Conversion(gv1)*CauchyTransform()
    r2 = Conversion(gv2)*CauchyTransform() ⊞ BoundaryValue(1,gv2)*CauchyTransform()
    Cp = r1 ⊘ r2

    M = Multiplication(g)

    r1 = M*BoundaryValue(-1,gv1)*CauchyTransform() ⊞ M*Conversion(gv1)*CauchyTransform()
    r2 = M*Conversion(gv2)*CauchyTransform() ⊞ M*BoundaryValue(-1,gv2)*CauchyTransform()
    MCm = r1 ⊘ r2
    S = Cp - MCm

    cS = S*(s1 ⊕ s2)
    u = \(cS,[f,f],100)
    truth = 1.0220013279200346
    approx = (CauchyTransform()*u[1])(1im) + (CauchyTransform()*u[2])(1im) + 1
    @test abs(truth - approx) < 1e-10

    g = x -> 1 + 0.5exp(-30x^4)
    f = x -> g(x) - 1
    gd1 = JacobiMappedInterval(-1,0,0.0,0.0)
    gd2 = JacobiMappedInterval(0,1im,0.0,0.0)
    s1 = Jacobi(0.0,0.0,gd1)
    s2 = Jacobi(0.0,0.0,gd2)

    dgd1 = DirectedLobattoMappedInterval(-1,0)
    dgd2 = DirectedLobattoMappedInterval(0,1im)
    gv1 = GridValues(dgd1)
    gv2 = GridValues(dgd2)

    r1 = BoundaryValue(1,gv1)*CauchyTransform() ⊞ Conversion(gv1)*CauchyTransform()
    r2 = Conversion(gv2)*CauchyTransform() ⊞ BoundaryValue(1,gv2)*CauchyTransform()
    Cp = r1 ⊘ r2
    M = Multiplication(g)

    r1 = M*BoundaryValue(-1,gv1)*CauchyTransform() ⊞ M*Conversion(gv1)*CauchyTransform()
    r2 = M*Conversion(gv2)*CauchyTransform() ⊞ M*BoundaryValue(-1,gv2)*CauchyTransform()
    MCm = r1 ⊘ r2
    S = Cp - MCm

    cS = S*(s1 ⊕ s2)
    u = \(cS,[f,f],100)

    truth = 0.9758397860028803 + 0.016084935690312958*im
    approx = (CauchyTransform()*u)(1.0) + 1
    @test abs(truth - approx) < 1e-10
end

