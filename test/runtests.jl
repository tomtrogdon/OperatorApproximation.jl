using OperatorApproximation, SpecialFunctions, LinearAlgebra, Memoize
using Test

@testset "OperatorApproximation.jl: Basic tests" begin
    gd = UltraMappedInterval(-10.0,10.0,0.0)
    sp = Ultraspherical(0,gd)
    gv = GridValues(gd)
    f = x -> exp(-x^2)
    fsp = BasisExpansion(f,sp,150)
    @test abs(fsp(.3) - f(.3)) < 1e-10

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

@testset "OperatorApproximation.jl: Hermite testing" begin
    gd = HermiteRealAxis()
    sp1 = HermiteFun(gd)
    sp2 = Erf(gd)
    D = Derivative()
    Op = (D ⊞ D)*(sp2 ⊕ sp1)
    f = x -> exp(-x^2/2)
    F = x -> (2pi)^(0.25)*exp(-x^2/2 + x^2/4)
    ff = CoefConversion(sp1)*BasisExpansion(F, HermitePoly(gd), 100)
    sol = \(Op,[ff],200)
    @test abs(sol(0.0) - sqrt(pi/2)) < 1e-12
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

@testset "OperatorApproximation.jl: Schrödinger scattering" begin
    @memoize function Sc(k,V)
        L = 20.0
        gdL = UltraMappedInterval(-L,0,0.0)
        gdR = UltraMappedInterval(0,L,0.0)
        spL = Ultraspherical(0.0,gdL)
        spR = Ultraspherical(0.0,gdR)
        D = Derivative()
        M = Multiplication(V)
        EL = Conversion(gdL |> GridValues)
        ER = Conversion(gdR |> GridValues)
        OpL = EL*D^2 - (2im*k)*(EL*D) + M*EL
        OmL = EL*D^2 + (2im*k)*(EL*D) + M*EL
        OpR = ER*D^2 - (2im*k)*(ER*D) + M*ER
        OmR = ER*D^2 + (2im*k)*(ER*D) + M*ER
        BL = Conversion(FixedGridValues([-L],gdL))
        BR = Conversion(FixedGridValues([L],gdR))
    
        OpL = (BL ⊘ (BL*D) ⊘ OpL)*spL
        JpL = OpL\([[0.0], [0.0], x -> -V(x)])
        dJpL = D*JpL
        OmL = (BL ⊘ (BL*D) ⊘ OmL)*spL
        JmL = OmL\([[0.0], [0.0], x -> -V(x)])
        dJmL = D*JmL
    
        OpR = (BR ⊘ (BR*D) ⊘ OpR)*spR
        JpR = OpR\([[0.0], [0.0], x -> -V(x)])
        dJpR = D*JpR
        OmR = (BR ⊘ (BR*D) ⊘ OmR)*spR
        JmR = OmR\([[0.0], [0.0], x -> -V(x)])
        dJmR = D*JmR
    
        A = [JpL(0.0)+1  JmL(0.0)+1;
            dJpL(0.0)-im*k*(JpL(0.0)+1) dJmL(0.0)+im*k*(JmL(0.0)+1)]
        B = [JpR(0.0)+1  JmR(0.0)+1;
            dJpR(0.0)-im*k*(JpR(0.0)+1) dJmR(0.0)+im*k*(JmR(0.0)+1)]
        B\A
    end
    ρ(k,V) = Sc(k,V)[2,1]/Sc(k,V)[1,1]
    # Without poles
    V = x -> -exp(-x^2)
    p = x -> ( abs(x) < 1e-12 ? -1.0 + 0.0im : ρ(x,V) )
    p̄ = x -> -conj(p(conj(x)))
    τ = x -> 1 + p(x)*p̄(x)
    JJ = [τ p̄; p x->1]
    Γ = [-10.0 0; 0 10.0]
    J = [JJ, JJ]
    rhp = RHP(Γ,J)
    rhsolver = RHSolver(rhp);
    sol = rhsolver([1 1], 300)
    U = CauchyTransform()*sol
    β = moment(sol[1],1)*1im/(2pi)
    α = moment(sol[1],0)*1im/(2pi)
    @test abs(-2(2β - α^2) - V(0.0)) < 1e-9
    
    # With pole
    V = x -> 1.3exp(-x^2)
    L = 23
    gd = PeriodicMappedInterval(-L,L)
    sp = Fourier(gd)
    D = Derivative()
    M = Multiplication(x -> V(x))
    Λ = eigen((-D^2 - M)*sp,1000)
    k = sqrt(Λ.values[1] |> real |> complex)
    ϵ = 1e-12
    h = x -> Sc(1im*x,V)[1,1]
    c = imag(h(imag(k) + 1im*ϵ))/ϵ # this only works because the data is decreasing so rapidly
    c = 1im/c
    
    p = x -> ( abs(x) < 1e-12 ? -1.0 + 0.0im : ρ(x,V) )
    p̄ = x -> -conj(p(conj(x)))
    τ = x -> 1 + p(x)*p̄(x)
    JJ = [τ p̄; p x->1]
    Γ = [-10.0 0; 0 10.0]
    J = [JJ, JJ]
    poles = [k, -k]
    Rp = [0 0; c 0]
    Rm = Rp' |> Matrix
    R = [Rp, Rm]
    rhp = RHP(Γ,J,poles,R)
    rhsolver = RHSolver(rhp);
    sol = rhsolver([1 1], 300)
    β = moment(sol[1],1)*1im/(2pi)
    α = moment(sol[1],0)*1im/(2pi)
    @test abs(-2(2β - α^2) - V(0.0)) < 1e-9
end

@testset "OperatorApproximation.jl: OscRational testing" begin
    #Function interpolation
    N = 200;
    α = 0;
    gd = RationalRealAxis()
    sp = OscRational(gd,α);
    f = x -> exp(-x^2)
    ff = BasisExpansion(f,sp,N)
    x_test = 0.1
    @test abs(f(x_test) - ff(x_test)) < 1e-10

    #Multplication
    M = Multiplication(ff)
    g = M*M*ff
    @test abs(g(x_test)-ff(x_test)^3) < 1e-10

    #Derivative
    D = Derivative()
    h = D*ff
    df = x-> -2*x*exp(-x^2)
    @test abs(df(x_test) - h(x_test)) < 1e-10

    #Alternative way to test derivative
    Op  = D*ff.basis
    dc = Matrix(Op,N+2,N)*ff.c
    test = BasisExpansion(ff.basis, dc)
    @test abs(test(x_test)-df(x_test)) < 1e-10

    #Oscillatory tests

    #Function interpolation
    N = 201;
    α = 2;
    gd = RationalRealAxis()
    sp = OscRational(gd,α);
    f = x -> exp(-x^2)
    ff = BasisExpansion(f,sp,N)
    x_test = 0.1
    @test abs(f(x_test)*exp(1im*α*x_test) - ff(x_test)) < 1e-10

    #Multplication
    M = Multiplication(ff)
    g = M*M*ff
    @test abs(g(x_test)-ff(x_test)^3) < 1e-10

    #Derivative
    D = Derivative()
    h = D*ff
    df = x-> -2*x*exp(-x^2)*exp(1im*α*x) + exp(-x^2)*1im*α*exp(1im*α*x)
    @test abs(df(x_test) - h(x_test)) < 1e-10

    #Cauchy tests

    function CPquad(f,a,z)
        R = 10
        fff = x -> f(x)*exp(1im*x*a)/(2im*pi)
        F = x -> f(x)*exp(1im*x*a)/(2im*pi)*1/(x-z)
        endpts = [-R, z-1, z-1im, z+1, R]
        s = 0
        for i = 1:4
            gd = JacobiMappedInterval(endpts[i],endpts[i+1],0,0)
            sp = Jacobi(0,0,gd)
            fff = BasisExpansion(F,sp)
            s = s + sum(fff)
        end
        s
    end

    function CMquad(f,a,z)
        R = 10
        fff = x -> f(x)*exp(1im*x*a)/(2im*pi)
        F = x -> f(x)*exp(1im*x*a)/(2im*pi)*1/(x-z)
        endpts = [-R, z-1, z+1im, z+1, R]
        s = 0
        for i = 1:4
            gd = JacobiMappedInterval(endpts[i],endpts[i+1],0,0)
            sp = Jacobi(0,0,gd)
            fff = BasisExpansion(F,sp)
            s = s + sum(fff)
        end
        s
    end

    #Nonoscillatory test
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

    

    #Oscillatory test (α > 0)
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


    #Oscillatory test (α > 0)
    α = -2;
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


    gd = RationalRealAxis();
    sp = OscRational(gd,0.0)
    ff = BasisExpansion(x -> exp(-x^2), sp, 300)
    spα = OscRational(gd,1.0)
    ffα = BasisExpansion(x -> exp(-x^2), spα, 300)
    fadd = ff + ffα;
    Cop = CauchyOperator(1)*fadd.basis;
    out = Cop*fadd

    @test abs(CPquad(x -> exp(-x^2),1,3) + CPquad(x -> exp(-x^2),0,3) - out(3)) < 1e-10

    #Integral test
    α = -2;
    gd = RationalRealAxis()
    sp = OscRational(gd,α);
    f = x -> exp(-x^2);
    ff = BasisExpansion(f,sp)
    int_test = Base.sum(ff)
    int_true = sqrt(π)*exp(-(abs(α)^2)/4)
    @test abs(int_test-int_true) < 1e-10

    #Testing Multiplication with two different α values for basis
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
    @test abs(h(x_test)-ff(x_test)*gg(x_test)) < 1e-10

    #Testing conjugate of a basis expansion
    α = 0;
    gd = RationalRealAxis()
    sp = OscRational(gd,α);
    f = x -> x/(x+1im)
    ff = BasisExpansion(f,sp)
    x_test = 5.67
    ff2 = Base.conj(ff)
    conj_f = x -> x/(x-1im)
    @test abs(ff2(x_test) - conj_f(x_test)) < 1e-10

    #Testing dot product of BasisExpansion{OscRational}
    α = 2;
    gd = RationalRealAxis()
    sp = OscRational(gd,α)
    f = x -> exp(-x^2)
    ff = BasisExpansion(f,sp)

    α = -3;
    sp = OscRational(gd,α)
    g = x -> exp(-2x^2)
    gg = BasisExpansion(g,sp)

    @test abs(dot(ff,gg) - sqrt(pi/3)*exp(-25/12)) < 1e-10

    #Testing dot product of BasisExpansion{DirectSum}
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
    
    @test abs(dot(h1,h2) - 2*sqrt(pi/3)*exp(-25/12)) < 1e-10

    #Test sumdot version of dot product of BasisExpansion{DirectSum}
    α = 2;
    gd = RationalRealAxis()
    sp = OscRational(gd,α)
    f = x -> exp(-x^2)
    ff = BasisExpansion(f,sp)

    α2 = -3;
    sp2 = OscRational(gd,α2)
    g = x -> exp(-2x^2)
    gg = BasisExpansion(g,sp2)

    h1 = ff ⊕ gg
    h2 = gg ⊕ ff

    true_val = 2*sqrt(pi/3)*exp(-25/12) + sqrt(pi/2) + sqrt(pi)/2
    @test abs(sumdot(h1,h2) - true_val) < 1e-10

    #Test vector version of sumdot
    v1 = [h1,h2]
    v2 = [h1,h2]
    # Need to fix true_val
    #@test abs(sumdot(v1,v2)-true_val*4) < 1e-10
end

@testset "OperatorApproximation.jl: erfc" begin
    gd = RationalMappedAxis(3.0,0.0,pi/2)
    sp = OscRational(gd,0.0)
    f = BasisExpansion(x -> 2exp(x^2),sp)
    Cf = CauchyTransform()*f
    merfc = z -> (real(z) > 0 ? - exp(-z^2)*Cf(z) : 2 - exp(-z^2)*Cf(z))
    @test abs(erfc(.3) - merfc(.3)) < 1e-12
    @test abs(erfc(-.3) - merfc(-.3)) < 1e-12
end