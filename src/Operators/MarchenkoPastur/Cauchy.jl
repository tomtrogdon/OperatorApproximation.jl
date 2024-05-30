function *(C::CauchyTransform,domain::MarchenkoPastur)
    a = domain.GD.D.a
    b = domain.GD.D.b
    d = domain.d
    gd = MarchenkoPasturMappedInterval(a,b,d)
    range = Hardy(Exterior(gd)) # Just to ensure the right weight is encoded
    ConcreteOperator(domain,range,BasicBandedOperator{ℕ₊,ℕ₊}(0,0, (i,j) -> i == j && i >= 0 ? complex(1.0) : 0.0im ))
end