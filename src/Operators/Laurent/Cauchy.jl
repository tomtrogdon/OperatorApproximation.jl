function *(C::Cauchy{T},domain::Laurent) where T <: Int64
    if C.o == 1
        range = Hardy(Interior(domain.GD))
        ConcreteLazyOperator(domain,range,BasicBandedOperator{ℤ,ℕ₊}(0,0, (i,j) -> i == j && i >= 0 ? complex(1.0) : 0.0im ))
    elseif C.o == -1
        range = Hardy(Exterior(domain.GD))
        ConcreteLazyOperator(domain,range,BasicBandedOperator{ℤ,ℕ₋}(0,0, (i,j) -> i == j && i < 0 ? complex(-1.0) : 0.0im ))
    end
end