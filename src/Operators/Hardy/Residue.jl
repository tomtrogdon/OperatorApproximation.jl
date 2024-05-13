function *(B::Residue,b1::Hardy{T,S}) where {T <: Exterior, S <: Interval}
    Op = ZeroOperator{ℕ₊,ℕ₊}(0,0,x -> 0)
    ConcreteOperator(b1,B.range,Op)
end

function *(B::Residue,b1::Hardy{T,S}) where {T <: Exterior, S <: DiscreteDomain}
    if !(typeof(B.range) <: FixedGridValues)
        @error "Wrong range for residue operator"
        return
    end
    ran_pts = B.range.GD.grid
    dom_pts = b1.GD.grid
    n = length(dom_pts)
    m = length(ran_pts)
    A = repeat(reshape(dom_pts,1,:),m,1) - repeat(reshape(ran_pts,:,1),1,n)
    A = map( x -> abs(x) < 1e-14 ? 1.0 : 0.0, A) |> sparse |> dropzeros   
    ConcreteOperator(b1,B.range,FixedMatrix{ℕ₊,ℕ₊}(A))
end