module HelperFunctions
export available_landau_levels, process_eos_data, bound_kz_tilde, neutron_density, proton_density, electron_density
include("Constants.jl")
using QuadGK


"""
    available_landau_levels(spin, k_F, mass, B, T)
    
Calculate the Landau levels available to a particle, including spin and thermal effects.
"""
function available_landau_levels(spin, k_F, mass, B, T)
    k_max_squared = k_F^2 + (TEMP_STEPS*T)^2 + 2*TEMP_STEPS*T*sqrt(k_F^2+mass^2) + spin*PROTON_G_MINUS_2*ELEM_CHARGE*B
    floor(k_max_squared / (2*ELEM_CHARGE*B))
end


"""
    process_eos_data(data, mageos)
    
Take raw EOS data from eosmlk or mageosmlk and process it into a uniform format.
"""
function process_eos_data(data, mageos)
    if mageos
        proton_bare_energy, neutron_bare_energy, _, _, mu_e, _, _, M_p_star, M_n_star = data
        k_Fn = sqrt(neutron_bare_energy^2 - M_n_star^2)
        k_Fp = sqrt(proton_bare_energy^2 - M_p_star^2)
    else
        k_Fp, k_Fn, mu_e, _, _, M_p_star, M_n_star = data
    end
    
    (k_Fn, k_Fp, mu_e, M_n_star, M_p_star)
end


"""
    bound_kz_tilde(k_F, mass, n, spin, B, T)

Calculate the bounds on the phasespace integral.
"""
function bound_kz_tilde(k_F, mass, n, spin, B, T)
    bare_energy = sqrt(k_F^2 + mass^2)

    kz_lower_squared = (bare_energy - TEMP_STEPS*T)^2 - mass^2 - 2*n*ELEM_CHARGE*B + spin*PROTON_G_MINUS_2*ELEM_CHARGE*B
    kz_upper_squared = (bare_energy + TEMP_STEPS*T)^2 - mass^2 - 2*n*ELEM_CHARGE*B + spin*PROTON_G_MINUS_2*ELEM_CHARGE*B
    if kz_upper_squared < 0
        return (0, 0.0001)
    end
    if kz_lower_squared < 0
        return (0, sqrt(kz_upper_squared)/T)
    end
    (sqrt(kz_lower_squared)/T, sqrt(kz_upper_squared)/T)
end


"""
    fermi_dirac(E, µ, T)

Calculate the Fermi-Dirac distribution for a given ``E, µ, T``.
"""
function fermi_dirac(E, µ, T)
    if (E-µ)/T > 50 return 0 end;
    if (E-µ)/T < -50 return 1 end;
    return 1 / (1 + exp((E - µ) / T))
end

Ee(k_ze, n_e, B) = sqrt(ELECTRON_MASS^2 + k_ze^2 + 2 * n_e * ELEM_CHARGE * (B * GAUSS_TO_MEV2))
En(k_n, spin, B) = sqrt(k_n^2 + NEUTRON_MASS^2) + spin * NEUTRON_AMM * (B * GAUSS_TO_MEV2)
Ep(k_zp, n_p, spin, B) = sqrt(PROTON_MASS^2 + k_zp^2 + 2 * n_p * ELEM_CHARGE * (B * GAUSS_TO_MEV2) - spin * PROTON_G_MINUS_2 * ELEM_CHARGE * (B * GAUSS_TO_MEV2))


"""
    electron_density(µ_e, B, T)

Calculate the electron number density as a function of chemical potential.
"""
function electron_density(µ_e, B, T)
    µe_gtr_me = µ_e > ELECTRON_MASS ? 1 : 0
    n_max = round(Int, floor(((µe_gtr_me*µ_e + TEMP_STEPS*T)^2 - ELECTRON_MASS^2) / (2 * ELEM_CHARGE * (B * GAUSS_TO_MEV2))))
    result = 0
    for n_e in 0:n_max+1
        k_z_max = sqrt(max(0, (µe_gtr_me*µ_e + TEMP_STEPS*T)^2 - ELECTRON_MASS^2 - 2 * n_e * ELEM_CHARGE * (B * GAUSS_TO_MEV2)))
        level_contrib = ELEM_CHARGE * (B * GAUSS_TO_MEV2) / (4*pi^2) * quadgk(k_ze -> fermi_dirac(Ee(k_ze, n_e, B), µ_e, T), -k_z_max, k_z_max)[1]
        if n_e > 0
            level_contrib *= 2
        end
        result += level_contrib
    end
    return result
end


"""
    proton_density(µ_p, B, T)

Calculate the proton number density as a function of chemical potential.
"""
function proton_density(µ_p, B, T)
    if µ_p > PROTON_MASS
        n_max = round(Int, floor(((µ_p + TEMP_STEPS*T)^2 - PROTON_MASS^2) / (2 * ELEM_CHARGE * (B * GAUSS_TO_MEV2))))
    else
        n_max = round(Int, floor((2 * PROTON_MASS * TEMP_STEPS * T) / (2 * ELEM_CHARGE * (B * GAUSS_TO_MEV2))))
    end
    
    result = 0
    for n_p in 0:n_max+2
        if µ_p > PROTON_MASS
            k_z_max_up_sq = (µ_p + TEMP_STEPS*T)^2 - PROTON_MASS^2 - 2*n_p*ELEM_CHARGE * (B * GAUSS_TO_MEV2) + PROTON_G_MINUS_2*ELEM_CHARGE * (B * GAUSS_TO_MEV2)/2
            k_z_max_down_sq = (µ_p + TEMP_STEPS*T)^2 - PROTON_MASS^2 - 2*n_p*ELEM_CHARGE * (B * GAUSS_TO_MEV2) - PROTON_G_MINUS_2*ELEM_CHARGE * (B * GAUSS_TO_MEV2)/2
        else
            k_z_max_up_sq = 2*PROTON_MASS*TEMP_STEPS*T - 2*n_p*ELEM_CHARGE * (B * GAUSS_TO_MEV2) + PROTON_G_MINUS_2*ELEM_CHARGE * (B * GAUSS_TO_MEV2)/2
            k_z_max_down_sq = 2*PROTON_MASS*TEMP_STEPS*T - 2*n_p*ELEM_CHARGE * (B * GAUSS_TO_MEV2) - PROTON_G_MINUS_2*ELEM_CHARGE * (B * GAUSS_TO_MEV2)/2
        end
        k_z_max_up = sqrt(max(0, k_z_max_up_sq))
        k_z_max_down = sqrt(max(0, k_z_max_down_sq))

        result += ELEM_CHARGE * (B * GAUSS_TO_MEV2) / (4*pi^2) * quadgk(k_zp -> fermi_dirac(Ep(k_zp, n_p, 0.5, B), µ_p, T), -k_z_max_up, k_z_max_up)[1]
        if n_p > 0
            result += ELEM_CHARGE * (B * GAUSS_TO_MEV2) / (4*pi^2) * quadgk(k_zp -> fermi_dirac(Ep(k_zp, n_p, -0.5, B), µ_p, T), -k_z_max_down, k_z_max_down)[1]
        end
    end
    return result
end


"""
    neutron_density(µ_n, B, T)

Calculate the neutron number density as a function of chemical potential.
"""
function neutron_density(µ_n, B, T)
    if µ_n > NEUTRON_MASS
        k_max_up_sq = (µ_n + TEMP_STEPS*T + NEUTRON_AMM * (B * GAUSS_TO_MEV2) / 2)^2 - NEUTRON_MASS^2
        k_max_down_sq = (µ_n + TEMP_STEPS*T - NEUTRON_AMM * (B * GAUSS_TO_MEV2) / 2)^2 - NEUTRON_MASS^2
    else
        k_max_up_sq = 2*NEUTRON_MASS * (TEMP_STEPS*T + NEUTRON_AMM * (B * GAUSS_TO_MEV2) / 2)
        k_max_down_sq = 2*NEUTRON_MASS * (TEMP_STEPS*T - NEUTRON_AMM * (B * GAUSS_TO_MEV2) / 2)
    end

    result = 0
    
    k_max_up = sqrt(max(0, k_max_up_sq))
    k_max_down = sqrt(max(0, k_max_down_sq))

    result += 1 / (2*pi^2) * quadgk(k_n -> k_n^2 * fermi_dirac(sqrt(k_n^2 + NEUTRON_MASS^2) + NEUTRON_AMM * (B * GAUSS_TO_MEV2) / 2, µ_n, T), 0, k_max_up)[1]
    result += 1 / (2*pi^2) * quadgk(k_n -> k_n^2 * fermi_dirac(sqrt(k_n^2 + NEUTRON_MASS^2) - NEUTRON_AMM * (B * GAUSS_TO_MEV2) / 2, µ_n, T), 0, k_max_down)[1]

    return result
end

end # module HelperFunctions