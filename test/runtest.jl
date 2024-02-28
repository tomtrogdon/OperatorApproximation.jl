using OperatorApproximation, SpecialFunctions
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
    @test abs(ff(3.0)-sin(3.0)) < 1e-14
end


@testset "OperatorApproximation.jl: Airy equation" begin
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
    setgrid(gv)
    u = [lbdry;rbdry;Op]\[[airyai(-R)];[airyai(R)]; x->0]
    @test abs(u(0) - airyai(0)) < 1e-10
end