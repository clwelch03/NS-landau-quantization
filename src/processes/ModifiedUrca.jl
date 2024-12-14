module ModifiedUrca
export modified_urca_rate
include("../utils/Constants.jl")


"""
    modified_urca_rate(k_Fn, k_Fp, m_n_star, m_p_star, T; cgs=true)

Calculate the modified Urca emissivity.
"""
function modified_urca_rate(k_Fn, k_Fp, m_n_star, m_p_star, T; cgs=true)
    m_n_landau = sqrt(k_Fn^2 + m_n_star^2)
    m_p_landau = sqrt(k_Fp^2 + m_p_star^2)
    n_p = k_Fp^3/(3*pi^2)
    n_0 = 1228857 # wow this sucks
    
    mass_ratios = (m_n_landau/m_n_star)^3 * (m_p_landau/m_p_star)
    density_ratio = (n_p/n_0)^(1/3)
    T_K9 = T * 11.61
    emissivity_cgs = 8.1*10^21 * mass_ratios * density_ratio * T_K9^8 * ALPHA_N * BETA_N
    if !cgs
        return emissivity_cgs * ERG_TO_MEV * RECIP_CM_TO_MEV^3 * RECIP_S_TO_MEV
    else
        return emissivity_cgs
    end
end

end # module MUrca