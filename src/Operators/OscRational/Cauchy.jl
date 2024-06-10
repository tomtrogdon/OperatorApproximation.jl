#inverse Mobius transform
Tm1 = z-> (1/1im)*((z+1)/(z-1))

#theta_vec
θ_vals = shift_mgrid(500,π)

#k_vec
k_vals = complex(zeros(length(θ_vals)-1))
for i=2:length(θ_vals)
    k_vals[i-1] = Tm1(exp(θ_vals[i]))
end

#define z
z = (k,σ)-> (-2im.*σ)./(k.+(σ.*1im))

#three term recurrence for generalized Lauguerre polynomials
function L_poly_x(x,α,nm1)
    num_polys = nm1+2
    L_poly = zeros(num_polys,1)
    if num_polys == 1
        L_poly[1] = 0
    elseif num_polys == 2
        L_poly[1] = 0
        L_poly[2] = 1
    else
        L_poly[1] = 0 #L_-1
        L_poly[2] = 1 #L_0
        L_poly[3] = 1+α-x #L_1
        L_iter = 3
        for n=1:num_polys-3
            L_np1α = ((((2*n)+1+α-x)*L_poly[L_iter])-((n+α)*L_poly[L_iter-1]))/(n+1) #L_n+1
            L_poly[n+3] = L_np1α
            L_iter += 1
        end
    end
    return L_poly
end

function M(k,σ,j)
    Mσ = complex(zeros(length(k),j))
    for h=1:length(k)
        Mσ[h,1] = 0 #r_0
        for i=1:j-1
            r_jα = ((1+z(k[h],σ))*Mσ[h,i]) + (z(k[h],σ)*(L_poly_x(2*abs(k[h]),1,i-1)[end] - L_poly_x(2*abs(k[h]),2,i-2)[end]))
            Mσ[h,i+1] = r_jα
        end
    end
    return Mσ
end

# #let's pseudocode some shit:
# function Cauchy(α,f,k_vals)
#     if α > 0
#         term1 = expansion of f in the Rational basis (How do I do this??)
#         term2 = R_{2m+1,0}*v(k_vals) (How exactly do I just extract R_{2m+1}??) (What is m, total number of coeffs???)
#         if want C+
#             return term1 + term2
#         else
#             return C+-I (Do I need to do anything special with this??)
#         end
#     else
#         if want C+
#             return -R_{2m+1,0}*v(k_vals)
#         else
#             return C+-I
#         end
#     end
# end
