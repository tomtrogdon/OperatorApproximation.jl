#This is creating the basis itself; it only depends on the interval that you want the grid
struct Rational <: Basis
    GD::GridAxis
end

#This defines the "domain" that the basis itself is defined on
cfd(sp::Rational) = ℤ #cfd=coefficient domain

#This defines the dimension of the basis
function dim(sp::Rational)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Rational
    nm = N₋(length(f.c))
    testconv(f.c[1:nm] |> reverse) && testconv(f.c[nm+1:end])
end

#Chops off coefficients less than 1e-15 and expands in Rational basis with those coeffs
function chop(f::BasisExpansion{T}) where T <: Rational
    nm = N₋(length(f.c)) #Int(floor(N/2))
    vm = chop(copy(f.c[1:nm]) |> reverse) |> reverse; #negative coeffs (chopping terms smaller than 1e-15)
    vp = chop(copy(f.c[nm+1:end])); #positive coeffs (chopping terms smaller than 1e-15)
    k = max(length(vm),length(vp)) #take max of length of vm or vp
    vm = pad(vm |> reverse,k) |> reverse #pad beginning with zeros to be equal to max length
    vp = pad(vp,k) #pad end with zeros to be equal to max length
    BasisExpansion(f.basis,vcat(vm,vp)) #expand in Rational basis using truncated coeffs with 1e-15 cutoff
end

#stealing code from Tom's 586 FFT tutorial
mfftshift = x -> circshift(fftshift(x), isodd(length(x)) ? 1 : 0)
mfft = x -> fftshift(fft(fftshift(x),1)) # fft(x,1) is used so that
# when we operate on matrices below, the behavior is as desired.
mifft = x -> mfftshift(ifft(mfftshift(x),1))
mgrid = (n,L) -> -L .+ 2*L*(0:n-1)/n
shift_mgrid = (n,L) -> (((-L .+ 2*L*(0:n-1)/n)./(2*L)).+(1/2)).*(2*π)
Pgrid = n -> 2*(0:n-1)/n .- 1

N₋ = N -> convert(Int64,floor(N/2))
N₊ = N -> convert(Int64,floor((N-1)/2))

#discrete Fourier transform
function mdft(v)
    if iseven(length(v))
        return mfft(v)/length(v)
    else # all this work for odd...
        m = length(v)
        mm = m ÷ 2
        σ = 1im*pi/m
        rot = exp(mm*σ)
        σ = exp(σ)
        w = mfft(v)
        @inbounds for i = 1:m
            w[i] *= rot
            rot /= σ
        end
        return w/m
    end
end

#inverse discrete Fourier transform
function midft(v)
    if iseven(length(v))
        return mifft(v*length(v))
    else
        w = copy(v)
        m = length(v)
        mm = m ÷ 2
        σ = 1im*pi/m
        rot = exp(mm*σ)
        σ = exp(σ)
        @inbounds for i = 1:m
            w[i] /= rot
            rot /= σ
        end
        return mifft(w*m)
    end
end

function (P::BasisExpansion{Rational})(k::Number,α::Number) # Horner's method
    x = P.basis.GD.D.imap(k)
    x = ((x.-1im)./(x.+1im))
    m = length(P.c)
    mm = convert(Int64,floor( m/2 )) #basically N₋(N);
    y = x^(-mm); 
    sum = P.c[1]*y
    for i = 2:length(P.c) # probably should do Horner from the zero mode outwards
        y  =  y.*x
        sum += P.c[i]*y
    end
    return exp.(1im*k*α).*sum
end
