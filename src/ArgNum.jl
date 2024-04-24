struct ArgNum <: Number
    z::ComplexF64
    ρ::Float64
    θ::Float64
    function ArgNum(z::Number,ρ::Float64,θ::Float64)
        new(z |> complex,ρ,θ |> mmod)
    end
end

function real(z::ArgNum)
    real(z.z)
end

function complex(z::ArgNum)
    z.z
end

function imag(z::ArgNum)
    imag(z.z)
end

function abs(z::ArgNum)
    abs(z.z)
end

function *(a::Number,z::ArgNum)
    ArgNum(a*z.z,z.ρ*abs(a),z.θ + angle(a))
end

function *(z::ArgNum,a::Number)
    a*z
end

function -(z::ArgNum)
    (-1)*z
end

function /(a::Number,z::ArgNum)
    ArgNum(a/z.z,abs(a)/z.ρ,-z.θ + angle(a))
end
function /(z::ArgNum,a::Number)
    ArgNum(z.z/a,z.ρ/abs(a),z.θ - angle(a))
end

function +(a::Number,z::ArgNum)
    ArgNum(z.z + a,z.ρ,z.θ)
end

function +(z::ArgNum,a::Number)
    ArgNum(z.z + a,z.ρ,z.θ)
end

function -(a::Number,z::ArgNum)
    a + (-z)
end

function -(z::ArgNum,a::Number)
    z + (-a)
end

function mmod(θ::Float64)
    angle(exp(1im*θ))
end

function log(z::ArgNum)
    if abs(z.z) < 1e-14
        return log(z.ρ) + 1im*mmod(z.θ) # return "finite part"
    # elseif imag(z.z) ≈ 0.0 && real(z.z) < 0
    #     return log(abs.(z.z)) + 1im*π*sign(mmod(z.θ))
    else
        return log(z.z)
    end
end