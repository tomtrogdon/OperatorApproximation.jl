function *(S::Shift,domain::Laurent)
    range = domain
    lb = S.order >= 0 ? abs.(S.order) : 0
    ub = S.order <= 0 ? abs(S.order) : 0
    ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(lb,ub, (i,j) -> i == j + S.order ? 1.0 : 0im ))
end