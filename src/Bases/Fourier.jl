struct Fourier <: Basis
    GD::GridInterval
end
####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
function dim(sp::Fourier)
    Inf
end

function pad(f::BasisExpansion{T},N) where T <: Fourier
    nm = N₋(length(f.c))
    new_nm = N₋(N)
    vm = pad(copy(f.c[1:nm]) |> reverse, new_nm) |> reverse;
    vp = pad(copy(f.c[nm+1:end]), N - new_nm);
    BasisExpansion(f.basis,vcat(vm,vp))
end

function testconv(f::BasisExpansion{T}) where T <: Fourier
    nm = N₋(length(f.c))
    testconv(f.c[1:nm] |> reverse) && testconv(f.c[nm+1:end])
end

function chop(f::BasisExpansion{T}) where T <: Fourier
    nm = N₋(length(f.c))
    vm = chop(copy(f.c[1:nm]) |> reverse) |> reverse;
    vp = chop(copy(f.c[nm+1:end]));
    k = max(length(vm),length(vp))
    vm = pad(vm |> reverse,k) |> reverse
    vp = pad(vp,k)
    BasisExpansion(f.basis,vcat(vm,vp))
end
####################################
####################################
####################################


mfftshift = x -> circshift(fftshift(x), isodd(length(x)) ? 1 : 0)
mfft = x -> fftshift(fft(fftshift(x),1)) # fft(x,1) is used so that
# when we operate on matrices below, the behavior is as desired.
mifft = x -> mfftshift(ifft(mfftshift(x),1))
Pgrid = n -> 2*(0:n-1)/n .- 1

N₋ = N -> convert(Int64,floor(N/2))
N₊ = N -> convert(Int64,floor((N-1)/2))

function mdft(v)
    if iseven(length(v))
        return mfft(v)/length(v)
    else # all this work for odd...
        m = length(v)
        mm = m ÷ 2
        σ = 1im*pi/m
        rot = exp(mm*σ)
        display(angle(rot))
        σ = exp(σ)
        w = mfft(v)
        @inbounds for i = 1:m
            w[i] *= rot
            rot /= σ
        end
        return w/m
    end
end

function midft(v)
    if iseven(length(v))
        return mifft(v*length(v))
    else
        w = copy(v)
        m = length(v)
        mm = m ÷ 2
        σ = 1im*pi/m
        rot = exp(mm*σ)
        display(angle(rot))
        σ = exp(σ)
        @inbounds for i = 1:m
            w[i] /= rot
            rot /= σ
        end
        return mifft(w*m)
    end
end

function (P::BasisExpansion{Fourier})(X::Number) # Horner's method
    x = P.basis.GD.D.imap(X)
    m = length(P.c)
    mm = convert(Int64,floor( m/2 ))
    ex = exp.(-1im*pi*mm*x)
    ex1 = exp.(1im*pi*x)
    sum = P.c[1]*ex
    for i = 2:length(P.c) # probably should do Horner from the zero mode outwards
        ex  =  ex.*ex1
        sum += P.c[i]*ex
    end
    return sum
end

