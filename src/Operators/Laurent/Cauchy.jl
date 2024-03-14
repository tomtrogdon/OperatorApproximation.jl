function *(C::Cauchy{T},domain::Laurent) where T <: Int64
    if C.o == 1
        range = domain
        dom = domain.GD.D
        ConcreteLazyOperator(domain,range,BasicBandedOperator(BI,0,0, (i,j) -> i == j && i >= 0 ? complex(1.0) : 0.0im ))
    elseif C.o == -1
        range = domain
        dom = domain.GD.D
        ConcreteLazyOperator(domain,range,BasicBandedOperator(BI,0,0, (i,j) -> i == j && i < 0 ? complex(-1.0) : 0im ))
    end
end