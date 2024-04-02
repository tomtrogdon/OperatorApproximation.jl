function RHdomain(endpoints::Matrix)
     temp = [(endpoints[i,1],endpoints[i,2]) for i=1:size(endpoints,2)]
     return ⊕([Legendre(t...) for t in temp]...)
end

function _rhrange(D::Basis)
    a = D.GD.D.a
    b = D.GD.D.b
    DirectedLobattoMappedInterval(a,b)
end

function RHrange(D::DirectSum)
    ⊕(GridValues.(_rhrange.(D.bases))...)
end

function RHmult(Js::Function) 

end

function RHmult(Js::Vector{T}) where T # J is a vector of scalar-valued functions

end

function RHmult(J::Vector{T}) where T <: Matrix # J is a vector of matrices of scalar-valued functions
    m = size(J[1],1)
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            g = [ JJ[i,j] for JJ in J]
            Js[j,i] = BlockDiagonalAbstractOperator(Multiplication.(g))
        end
    end
    Js
end

function RHrhs(J::Vector{T},c) where T <: Matrix # J is a vector of matrices of scalar-valued functions
    m = size(J[1],1)
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            Js[j,i] = [ z -> c*(ComplexF64.(map(x -> x(z),JJ)) - I) for JJ in J]
        end
    end
    Js
end