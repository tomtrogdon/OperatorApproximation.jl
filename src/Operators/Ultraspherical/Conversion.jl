function isconvertible(b1::Ultraspherical,b2::DiscreteBasis)
    iscompatible(b1.GD,b2.GD)
end

function isconvertible(b1::Ultraspherical,b2::Ultraspherical)
    iscompatible(b1.GD,b2.GD) && b2.λ >= b1.λ && mod(b1.λ - b2.λ,1) ≈ 0
end

function conversion(b1::Ultraspherical,b2::GridValues)
    basegrid =  n -> b2.GD.grid(n)
    # In principle, we need to do this:
    # gridfun = n -> b1.GD.D.imap(b2.GD.D.map(basegrid(n)))
    # but we are checking that the two grid domains are compatible
    # and currently this forces the composition of the maps to
    # be the identity
    a, b = Jacobi_ab(b1.λ - 1/2, b1.λ - 1/2)
    Op = OPEvaluationOperator(basegrid,a,b)
    ConcreteOperator(b1,b2,Op)
end

function conversion(b1::Ultraspherical,b2::FixedGridValues)
    # See conversion remark above.
    a, b = Jacobi_ab(b1.λ - 1/2, b1.λ - 1/2)
    Op = FixedGridOPEvaluationOperator(b2.pts,a,b)
    ConcreteOperator(b1,b2,Op)
end

function _s(j,λ)
    if j == 0
        return 1.0
    end
    0.5*((λ + 1)/(2*λ +1))*((j + 2λ)/(j + λ))*((j+2λ+1)/(j+λ+1)) |> sqrt
end

function _t(j,λ)
    0.5*((j-1)/(j + λ))*((j)/(j + λ - 1))*((λ+1)/(2λ+1)) |> sqrt
end

function _conv_ultra(λ)
    function out(i,j)
        if i == j
            return _s(j-1,λ)
        elseif j == i + 2
            return -_t(j-1,λ) 
        else
            return 0
        end
    end
end

function conversion(b1::Ultraspherical,b2::Ultraspherical)
    # λ0 = b1.λ
    # sp0 = Ultraspherical(λ0 + 1, b2.GD)
    # L0 = ConcreteOperator(b1,sp0,BasicBandedOperator(0,2,(i,j) -> _conv_ultra(i,j,λ0)))
    # λ0 += 1
    # while λ0 - b2.λ !≈ 0
    #     L0 = ConcreteOperator(b1,sp0,BasicBandedOperator(0,2,(i,j) -> _conv_ultra(i,j,λ0)))*L0
    #     λ0 += 1
    # end
    # L0
    if b1.λ ≈ b2.λ
         return ConcreteOperator(b1,b2,BasicBandedOperator{ℕ₊,ℕ₊}(0,0,(i,j) -> Float64(i == j)))
    end

    λ0 = b1.λ
    L0 = BasicBandedOperator{ℕ₊,ℕ₊}(0,2,_conv_ultra(λ0))
    λ0 += 1
    while !(λ0 - b2.λ ≈ 0)
        L0 = BasicBandedOperator{ℕ₊,ℕ₊}(0,2,_conv_ultra(λ0))*L0
        λ0 += 1
    end
    ConcreteOperator(b1,b2,L0)
end