function *(D::FourierTransform,domain::OscRational)
    o = D.o
    if abs.(domain.GD.D.θ) > 0
        @error "FourierTransform not supported on domain"
        return
    end
    α = o*domain.α # now only works for real axis
    if domain.GD.D == RealAxis()
        σ = 1.0
    else
        σ = 1/domain.GD.D.σ
    end
    dom_p = MappedSemiAxis(σ,α,0.0)
    gd_p = LaguerreSemiAxis(dom_p,1.0)
    sp_p = LaguerreFun(gd_p,1.0)

    dom_m = MappedSemiAxis(σ,α,pi)
    gd_m = LaguerreSemiAxis(dom_m,1.0)
    sp_m = LaguerreFun(gd_m,1.0)

    dom_o = Grid([α])
    gd_o = FixedGridValues(dom_o)
    range = sp_p ⊕ sp_m ⊕ gd_o




    # if D.order == 1
    #     range = domain
    #     ConcreteOperator(domain,range,BasicBandedOperator{ℤ,ℤ}(1,1,(i,j) -> DerivMatrix(i,j,α)))
    # else
    #     Derivative(D.order-1)*(Derivative(1)*domain)
    # end
end