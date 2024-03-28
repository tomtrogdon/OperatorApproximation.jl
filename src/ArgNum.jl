struct ArgNum <: Number
    z::ComplexF64
    θ::Float64
    function ArgNum(z::Number,θ::Float64)
        new(z |> complex,θ |> mmod)
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
    ArgNum(a*z.z,z.θ + angle(a))
end

function *(z::ArgNum,a::Number)
    ArgNum(a*z.z,z.θ + angle(a))
end

function -(z::ArgNum)
    (-1)*z
end

function /(a::Number,z::ArgNum)
    ArgNum(a/z.z,-z.θ + angle(a))
end
function /(z::ArgNum,a::Number)
    ArgNum(z.z/a,z.θ - angle(a))
end

function +(a::Number,z::ArgNum)
    ArgNum(z.z + a,z.θ)
end

function +(z::ArgNum,a::Number)
    ArgNum(z.z + a,z.θ)
end

function -(a::Number,z::ArgNum)
    a + (-z)
end

function -(z::ArgNum,a::Number)
    z + (-a)
end

function length(a::ArgNum)
    1
end

function iterate(a::ArgNum)
    (a,nothing)
end

function iterate(a::ArgNum,nothing)
end

function mmod(θ::Float64)
    angle(exp(1im*θ))
end

function log(z::ArgNum)
    if abs(z.z) ≈ 0.0
        return 1im*mmod(z.θ) # return "finite part"
    elseif imag(z.z) ≈ 0.0 && real(z.z) < 0
        return log(abs.(z.z)) + 1im*π*sign(mmod(z.θ))
    else
        return log(z.z)
    end
end

function sqrt(z::ArgNum)
    if imag(z.z) ≈ 0.0 && real(z.z) <= 0
        Z = sign(z.θ)*sqrt(abs(z.z))
        return ArgNum(Z, z.θ - angle(Z))
    else
        Z = sqrt(z.z)
        return ArgNum(Z,z.θ - angle(Z))
    end
end