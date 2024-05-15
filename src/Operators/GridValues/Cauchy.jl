function *(C::CauchyTransform,domain::FixedGridValues)
    gd = domain.GD
    range = Hardy(Exterior(gd)) # Just to ensure the right weight is encoded
    ConcreteOperator(domain,range,BasicBandedOperator{ℕ₊,ℕ₊}(0,0, (i,j) -> i == j && i >= 0 ? complex(1.0) : 0.0im ))
end

function poleres_cauchy(ps,z::Number)
    sv = (ps .- z)
    flag = abs.(sv) .< 1e-14
    sv[flag] .= ones(sum(flag))
    sv .= 1.0./sv
    sv[flag] .= zeros(sum(flag))
    1/(2im*pi)*sv
end

function poleres_cauchy(ps,z::Vector)  # vectorize!
    A = zeros(ComplexF64,length(z),length(ps))
    for i = 1:length(z)
        A[i,:] .= poleres_cauchy(ps,z[i])
    end
    A
end