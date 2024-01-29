struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector
end

function plot(f::BasisExpansion;dx = 0.01)
    x = -1:dx:1
    x = f.basis.GD.D.map.(x)
    plot(x,f.(x))
end

function *(Op::AbstractOperator,f::BasisExpansion)
    Opc = Op*f.basis
    Opc*f
end

function *(Op::ConcreteLazyOperator,f::BasisExpansion)
    BasisExpansion(Op.range,Op.L*f.c)
end