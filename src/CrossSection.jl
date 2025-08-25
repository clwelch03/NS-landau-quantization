module CrossSection

export cross_section, cross_section_zero_field

include("utils/HelperFunctions.jl")
include("utils/Constants.jl")
include("eos/EOS.jl")
include("utils/I_nr.jl")
 
using HCubature, Roots, SpecialFunctions, .EOS, .HelperFunctions, QuadGK, .I_nr

fermi_dirac_simple(x) = 1/(exp(x)+1) #TODO: all calls to this should eventually be phased out and turned into calls to fermi_dirac

"""
    m_red_low_density(proton_spin, neutron_spin, n_p, n_e, E_e, k_n, k_zp, k_ze, k_nu, cos_theta_nu, B)

Calculate the reduced matrix element at low density -- i.e. for nucleon capture.
"""
function m_red_low_density(proton_spin, neutron_spin, n_p, n_e, E_e, k_n, k_zp, k_ze, k_nu, cos_theta_nu, B)
    k_zn = k_zp + k_ze - k_nu*cos_theta_nu
    I_arg = (k_n^2 - k_zn^2) / (2*ELEM_CHARGE*B) # this is (w_perp^2 / 2eB)

    # see paper
    if proton_spin == 0.5 && neutron_spin == 0.5
        result = ( (G_V + G_A)^2 * (1 + k_ze/E_e) * (1 + cos_theta_nu) * I(n_e, n_p, I_arg)^2
                    + (G_V - G_A)^2 * (1 - k_ze/E_e) * (1 - cos_theta_nu) * I(n_e-1, n_p, I_arg)^2 )
    elseif proton_spin == 0.5 && neutron_spin == -0.5
        result = 4*G_A^2 * (1 + k_ze/E_e) * (1 - cos_theta_nu) * I(n_e, n_p, I_arg)^2
    elseif proton_spin == -0.5 && neutron_spin == 0.5
        result = 4*G_A^2 * (1 - k_ze/E_e) * (1 + cos_theta_nu) * I(n_e-1, n_p-1, I_arg)^2
    elseif proton_spin == -0.5 && neutron_spin == -0.5
        result = ( (G_V + G_A)^2 * (1 - k_ze/E_e) * (1 - cos_theta_nu) * I(n_e-1, n_p-1, I_arg)^2
                    + (G_V - G_A)^2 * (1 + k_ze/E_e) * (1 + cos_theta_nu) * I(n_e, n_p-1, I_arg)^2 )
    else
        @error "m_red_low_density: Either the proton spin $proton_spin or neutron spin $neutron_spin is not equal to ±0.5."
    end
    return result
end


"""
    cross_section_sign_sum(channel, neutron_spin, proton_spin, n_p, n_e, k_n, E_e, k_zp, k_ze, k_nu, cos_theta_nu, B)

Calculate the sign sum for the opacity calculation.
"""
function cross_section_sign_sum(channel, neutron_spin, proton_spin, n_p, n_e, k_n, E_e, k_zp, k_ze, k_nu, cos_theta_nu, B)
    sum_over_signs = 0
    k_znu = k_nu*cos_theta_nu

    # Step function based on whether we're doing n + ν -> p + e or p + \bar{ν} -> n + \bar{e}
    if channel=="n"
        for k_zp_sign in (-1, 1), k_ze_sign in (-1, 1)
            if k_n > abs(k_zp_sign*k_zp + k_ze_sign*k_ze - k_znu)
                sum_over_signs += m_red_low_density(proton_spin, neutron_spin, n_p, n_e, E_e, k_n, k_zp_sign*k_zp, k_ze_sign*k_ze, k_nu, cos_theta_nu, B)
            end
        end
    elseif channel=="p"
        for k_zp_sign in (-1, 1), k_ze_sign in (-1, 1)
            if k_n > abs(k_zp_sign*k_zp - k_ze_sign*k_ze + k_znu)
                sum_over_signs += m_red_low_density(proton_spin, neutron_spin, n_p, n_e, E_e, k_n, k_zp_sign*k_zp, k_ze_sign*k_ze, k_nu, cos_theta_nu, B)
            end
        end
    else
        @error "cross_section_sign_sum: channel not properly set as 'n' or 'p'"
    end
    sum_over_signs
end


"""
    _bound_kz_tilde_cross_section(n_p, n_e, k_nu, proton_spin, B, T, mu_p, mu_e)

Place bounds on kz_tilde for the cross-section phasespace integral.
"""
function _bound_kz_tilde_cross_section(n_p, n_e, k_nu, proton_spin, B, T, mu_p, mu_e)
    if mu_p < PROTON_MASS
        kzp_upper_squared = max(0, 2*PROTON_MASS*(k_nu + TEMP_STEPS*T + proton_spin*PROTON_G_MINUS_2*ELEM_CHARGE*B/PROTON_MASS - n_p*ELEM_CHARGE*B/PROTON_MASS))
    else
        kzp_upper_squared = max(0, 2*PROTON_MASS*(mu_p - PROTON_MASS + k_nu + TEMP_STEPS*T + proton_spin*PROTON_G_MINUS_2*ELEM_CHARGE*B/PROTON_MASS - n_p*ELEM_CHARGE*B/PROTON_MASS))
    end

    if mu_e < ELECTRON_MASS
        kze_upper_squared = max(0, (k_nu + TEMP_STEPS*T)^2 - (ELECTRON_MASS^2 + 2*n_e*ELEM_CHARGE*B))
    else
        kze_upper_squared = max(0, (mu_e + k_nu + TEMP_STEPS*T)^2 - (ELECTRON_MASS^2 + 2*n_e*ELEM_CHARGE*B))
    end

    (0, sqrt(kzp_upper_squared)/T), (0, sqrt(kze_upper_squared)/T)
end


"""
    cross_section_integrand(channel, k_zp_tilde, k_ze_tilde, n_p, n_e, mu_n, mu_p, mu_e, k_nu, cos_theta_nu, neutron_spin, proton_spin, B, T)

Compute the integrand for the cross-section calculation.
"""
function cross_section_integrand(channel, k_zp_tilde, k_ze_tilde, n_p, n_e, mu_n, mu_p, mu_e, k_nu, cos_theta_nu, neutron_spin, proton_spin, B, T, blocking)
    k_zp = k_zp_tilde * T
    k_ze = k_ze_tilde * T
    E_p = sqrt(PROTON_MASS^2 + k_zp^2 + 2*n_p*ELEM_CHARGE*B - PROTON_G_MINUS_2*ELEM_CHARGE*B*proton_spin)
    E_e = sqrt(ELECTRON_MASS^2 + k_ze^2 + 2*n_e*ELEM_CHARGE*B)
    E_n = channel == "n" ? E_p + E_e - k_nu : E_p + k_nu - E_e

    # Checks that the state is accessible
    if channel == "n"
        if E_e + E_p - k_nu >= maximum([NEUTRON_MASS, mu_n]) + TEMP_STEPS * T
            return 0
        end
    elseif channel == "p"
        if E_p + k_nu - E_e < NEUTRON_MASS + B * NEUTRON_AMM / 2
            return 0
        end
    end

    k_n_squared = (E_n + neutron_spin*NEUTRON_AMM*B)^2 - NEUTRON_MASS^2
    if k_n_squared < 0
        return 0
    else
        k_n = sqrt(k_n_squared)
    end
    
    # Fermi-Dirac factors, ignoring final state blocking for the electron
    if channel == "n"
        if blocking
            fermi_dirac_factors = fermi_dirac_simple((E_e+E_p-k_nu-mu_n)/T) * fermi_dirac_simple((-E_p+mu_p)/T) * fermi_dirac_simple((-E_e+mu_e)/T)
        else
            fermi_dirac_factors = fermi_dirac_simple((E_e+E_p-k_nu-mu_n)/T)
        end
    elseif channel == "p"
        if blocking
            fermi_dirac_factors = fermi_dirac_simple((E_e-E_p-k_nu+mu_n)/T) * fermi_dirac_simple((E_p-mu_p)/T) * fermi_dirac_simple((-E_e-mu_e)/T)
        else
            fermi_dirac_factors = fermi_dirac_simple((E_p-mu_p)/T)
        end
    else
        @error "cross_section_integrand: cross_section integrand channel not properly set as 'n' or 'p'"
    end
        
    sign_sum = cross_section_sign_sum(channel, neutron_spin, proton_spin, n_p, n_e, k_n, E_e, k_zp, k_ze, k_nu, cos_theta_nu, B)

    # Magic factor of 2 to match DQ
    2 * E_n*fermi_dirac_factors*sign_sum
end


"""
    cross_section_integral(channel, n_p, n_e, mu_n, mu_p, mu_e, k_nu, cos_theta_nu, neutron_spin, proton_spin, B, T)

Compute the integrand for the cross-section calculation.
"""
function cross_section_integral(channel, n_p, n_e, mu_n, mu_p, mu_e, k_nu, cos_theta_nu, neutron_spin, proton_spin, B, T, blocking, evals)
    kzp_tilde_bounds, kze_tilde_bounds = _bound_kz_tilde_cross_section(n_p, n_e, k_nu, proton_spin, B, T, mu_p, mu_e)

    integrand_func_wrapper(kz_tilde_vec) = cross_section_integrand(channel, kz_tilde_vec[1], kz_tilde_vec[2], n_p, n_e, mu_n,
                                                mu_p, mu_e, k_nu, cos_theta_nu, neutron_spin, proton_spin, B, T, blocking)
    
    integral_result, _ = hcubature(integrand_func_wrapper, (kzp_tilde_bounds[1], kze_tilde_bounds[1]), 
        (kzp_tilde_bounds[2], kze_tilde_bounds[2]); maxevals=evals) # huge number of maxevals so the integral actually works
    integral_result
end


"""
    cross_section(channel, n_B, Y_e, k_nu, cos_theta_nu, B, T)

Compute the integrand for the cross-section calculation.
"""
function cross_section(channel, n_B, Y_e, k_nu, cos_theta_nu, B, T; blocking = true, evals = 50000)
    B_MeV2 = B * GAUSS_TO_MEV2
    prefactor = (G_F^2 * COS2_CABIBBO_ANGLE * ELEM_CHARGE * B_MeV2 * T^2) / (64*pi^3)

    n_n = (1-Y_e) * n_B
    n_p = Y_e * n_B

    mu_n = find_zero(mu -> n_n - neutron_density(mu, B, T), (850, 1000))
    mu_p = find_zero(mu -> n_p - proton_density(mu, B, T), (850, 1000))
    mu_e = find_zero(mu -> n_p - electron_density(mu, B, T), (0, 200))
    println(mu_n, " ", mu_p, " ", mu_e)
    # Maximum Landau level attainable for proton (spin up/down) and electron
    # If spin splitting is disabled, n_max_p_up = n_max_p_down, and we can use
    # either one of them
    n_max_p_up = floor(PROTON_MASS/(ELEM_CHARGE*B_MeV2) * (k_nu + TEMP_STEPS*T + 0.25*PROTON_G_MINUS_2*ELEM_CHARGE*B_MeV2/PROTON_MASS - ELECTRON_MASS))
    n_max_p_down = floor(PROTON_MASS/(ELEM_CHARGE*B_MeV2) * (k_nu + TEMP_STEPS*T - 0.25*PROTON_G_MINUS_2*ELEM_CHARGE*B_MeV2/PROTON_MASS - ELECTRON_MASS))
    n_max_e = floor(((k_nu + TEMP_STEPS*T + 0.25*PROTON_G_MINUS_2*ELEM_CHARGE*B_MeV2/PROTON_MASS)^2 - ELECTRON_MASS^2) / (2*ELEM_CHARGE*B_MeV2))

    # Sum up the contributions from each pair of LLs
    landau_sum = 0
    for neutron_spin in (-0.5, 0.5), proton_spin in (-0.5, 0.5)
        for n_e in 0:n_max_e
            for n_p_up in 0:n_max_p_up
                landau_sum += cross_section_integral(channel, n_p_up, n_e, mu_n, mu_p, mu_e, k_nu, cos_theta_nu, neutron_spin, proton_spin, B_MeV2, T, blocking, evals)
            end
            for n_p_down in 0:n_max_p_down
                landau_sum += cross_section_integral(channel, n_p_down, n_e, mu_n, mu_p, mu_e, k_nu, cos_theta_nu, neutron_spin, proton_spin, B_MeV2, T, blocking, evals)
            end
        end
    end

    if channel == "n"
        prefactor * landau_sum * (197.3/10^13)^2/n_n
    elseif channel == "p"
        prefactor * landau_sum * (197.3/10^13)^2/n_p
    else
        @error "cross_section: channel not properly set as 'n' or 'p'"
    end
end


Δnp = 1.293 # mass gap


"""
    cross_section_approx(channel, k_nu, cos_theta_nu, B, T)

Compute the zero-field approximation, with first-order magnetic correction, from the Duan/Qian paper.
"""
function cross_section_approx(channel, k_nu, cos_theta_nu, B, T)
    np_sign = channel == "n" ? +1 : -1
    χ = channel == "n" ? NEUTRON_AMM * B * GAUSS_TO_MEV2/T : PROTON_AMM * B * GAUSS_TO_MEV2/T
    numerator_pm = channel == "n" ? G_V + G_A : G_V - G_A
    ϵ = χ * (2*numerator_pm*G_A)/(G_V^2+3*G_A^2) * cos_theta_nu
    E_e_0 = max(ELECTRON_MASS, k_nu + np_sign*Δnp)
    mu_νN_0 = G_F^2*COS2_CABIBBO_ANGLE/pi * E_e_0 * sqrt(E_e_0^2 - ELECTRON_MASS^2) * (197.3/10^13)^2
    mu_νN_1 = mu_νN_0 * (1 - (2*(G_V^2 - np_sign*2*(G_V+F_2)*G_A + 5*G_A^2)) / (G_V^2 + 3*G_A^2) * k_nu/NEUTRON_MASS)
    mu_νN_1 * (1+ϵ)
end

end # module sigma