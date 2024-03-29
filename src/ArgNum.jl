struct ArgNum
    z::ComplexF64
    θ::Float64
end

function *(a::Number,z::ArgNum)
    ArgNum(a*z.z,z.θ + angle(a))
end
function *(z::ArgNum,a::Number)
    ArgNum(a*z.z,z.θ + angle(a))
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
    ArgNum(z.z - a,z.θ)
end
function -(z::ArgNum,a::Number)
    ArgNum(z.z - a,z.θ)
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
        return 1im*mmod(z.θ)
    elseif imag(z.z) ≈ 0.0 && real(z.z) < 0
        return log(abs.(z.z)) + 1im*π*sign(mmod(z.θ))
    else
        return log(z.z)
    end
end

function sqrt(z::ArgNum)
    if imag(z.z) ≈ 0.0 && real(z.z) < 0
        return sign(mmod(z.θ))*sqrt(abs(z.z))
    else
        return sqrt(z.z)
    end
end