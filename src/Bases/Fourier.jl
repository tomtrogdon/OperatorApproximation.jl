struct Fourier <: Basis
    GD::GridDomain
end

# function BasisExpansion(f::Function,basis::Fourier,N::Integer)
#     Conversion(basis)*BasisExpansion(f,GridValues(basis.GD),N)
# end