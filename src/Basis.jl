abstract type Basis end

struct GridValues <: Basis
    GD::GridDomain
end

function iscompatible(GV1::GridValues,GV2::GridValues)
    GV1 == GV2
end