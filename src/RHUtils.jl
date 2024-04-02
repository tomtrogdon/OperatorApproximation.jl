function RHdomain(endpoints::Matrix)
     temp = [(endpoints[i,1],endpoints[i,2]) for i=1:size(endpoints,1)]
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
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Matrix{Any}(nothing,m,m)
    for i = 1:m
        for j = 1:m
            g = [ JJ[i,j] for JJ in J]
            Js[j,i] = BlockDiagonalAbstractOperator(Multiplication.(g))
        end
    end
    convert(Matrix{BlockDiagonalAbstractOperator},Js)
end

function RHrhs(J::Vector{T},c) where T <: Matrix # J is a vector of matrices of scalar-valued functions
    m = size(J[1],1) # m is size of RHP
    # length of J is # of contours
    Js = Vector{Any}(nothing,m)
    for i = 1:m
        Js[i] = [ z -> (c*(ComplexF64.(map(x -> x(z),JJ)) - I))[i] for JJ in J]
    end
    convert(Vector{Vector},Js)
end