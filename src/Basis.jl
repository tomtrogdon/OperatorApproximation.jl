abstract type Basis end

struct GridValues <: Basis
    GD::GridDomain
end
