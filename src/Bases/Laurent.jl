Lgrid = n -> exp.((2*(0:n-1)/n .- 1)*1im*pi )

struct Laurent <: Basis
    GD::GridCircle
end
####################################
#### REQUIRED TO BE IMPLEMENTED ####
####################################
cfd(sp::Laurent) = ℤ

function dim(sp::Laurent)
    Inf
end

function testconv(f::BasisExpansion{T}) where T <: Laurent
    nm = N₋(length(f.c))
    testconv(f.c[1:nm] |> reverse) && testconv(f.c[nm+1:end])
end

function chop(f::BasisExpansion{T}) where T <: Laurent
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

function (P::BasisExpansion{Laurent})(X::Number) # Horner's method
    x = P.basis.GD.D.imap(X)
    m = length(P.c)
    mm = convert(Int64,floor( m/2 ))
    ex = complex(float(x)).^(-mm)
    sum = P.c[1]*ex
    for i = 2:length(P.c) # probably should do Horner from the zero mode outwards
        ex  =  ex.*x
        sum += P.c[i]*ex
    end
    return sum
end