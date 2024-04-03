function *(C::CauchyTransform,b::BasisExpansion{T}) where T <: DirectSum
    ⊕([C*b[i] for i = 1:size(b,1)]...)
end

