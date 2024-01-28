struct BasisExpansion{T<:Basis}
    basis::T
    c::Vector
end

function plot(f::BasisExpansion)
    x = -1:0.01:1
    x = f.basis.GD.D.map.(x)
    plot(x,f.(x))
end