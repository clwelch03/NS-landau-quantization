module Urca
export direct_urca_rate, modified_urca_rate, density, qc_emissivity, find_µ_range_endpoints, rb, urca_rate_wrapper

include("eos/EOS.jl")
include("utils/HelperFunctions.jl")
include("utils/Constants.jl")
include("utils/I_nr.jl")

using HCubature, SpecialFunctions, .EOS, .HelperFunctions, PolyLog, Roots, .I_nr


"""
    m_red(spin, n_p, n_e, E_Fn, E_Fe, k_Fn, k_zp, k_ze, B, M_n_star, M_p_star)

Calculate the reduced matrix element based on state parameters.
"""
function m_red(proton_spin, n_p, n_e, E_Fn, E_Fe, k_Fn, k_zp, k_ze, B, M_n_star, M_p_star) # important: B must be in MeV^2!
    # In the no-spin-splitting case, we simply add the two parts of M_red.
    if proton_spin == 0
        return (m_red(0.5, n_p, n_e, E_Fn, E_Fe, k_Fn, k_zp, k_ze, B, M_n_star, M_p_star)
            + m_red(-0.5, n_p, n_e, E_Fn, E_Fe, k_Fn, k_zp, k_ze, B, M_n_star, M_p_star))
    end

    # Otherwise, we proceed with the calculation.
    M_L_p_star = sqrt(M_p_star^2 + k_zp^2 + 2*n_p*ELEM_CHARGE*B - ELEM_CHARGE*B*proton_spin*PROTON_G_MINUS_2)
    M_L_n_star = E_Fn

    k_zn = k_zp + k_ze # due to momentum conservation

    e_plus = 1 + k_ze/E_Fe
    e_minus = 1 - k_ze/E_Fe
    pz_plus = (1 + k_zp/(M_L_p_star + M_p_star))^2
    pz_minus = (1 - k_zp/(M_L_p_star + M_p_star))^2
    pB = (2*n_p*ELEM_CHARGE*B)/(M_L_p_star+M_p_star)^2

    I_arg = (k_Fn^2 - k_zn^2) / (2*ELEM_CHARGE*B) # this is (w_perp^2 / 2eB)

    I_ne_np = I(n_e, n_p, I_arg)
    I_neminus1_np = I(n_e-1, n_p, I_arg)
    I_ne_npminus1 = I(n_e, n_p-1, I_arg)
    I_neminus1_npminus1 = I(n_e-1, n_p-1, I_arg)

    # Differentiate between the proton spin-up and spin-down cases
    if proton_spin == 0.5
        gv_plus_ga_term = (16 * (G_V + G_A)^2 * M_L_n_star * (I_ne_np^2 * pz_minus * e_plus 
                                    + I_neminus1_npminus1^2 * pB * e_minus
                                    + 2 * I_ne_np * I_neminus1_npminus1 * (1-k_zp/(M_L_p_star+M_p_star))
                                        * (2*ELEM_CHARGE*B*sqrt(n_e*n_p))/((M_L_p_star+M_p_star)*E_Fe)))
    
        gv_minus_ga_term = (8 * (G_V - G_A)^2 * ((M_L_n_star - k_zn) * ((I_ne_np^2 * pz_plus * e_plus)
                                    + (I_ne_npminus1^2 * pB * e_plus))
                                    + (M_L_n_star + k_zn) * ((I_neminus1_np^2 * pz_plus * e_minus)
                                    + (I_neminus1_npminus1^2 * pB * e_minus))))
        
        gv2_minus_ga2_term = (16 * (G_V^2-G_A^2) * M_p_star * (I_ne_np^2 * (k_zp^2/(M_L_p_star + M_p_star)^2 - 1) * e_plus
                                    + I_neminus1_npminus1^2 * pB * e_minus
                                    - 2 * I_ne_np * I_neminus1_npminus1 * k_zp
                                        * (2*ELEM_CHARGE*B*sqrt(n_e*n_p))/((M_L_p_star+M_p_star)^2*E_Fe)))
    elseif proton_spin == -0.5
        gv_plus_ga_term = (16 * (G_V + G_A)^2 * M_L_n_star * (I_ne_np^2 * pB * e_plus 
                                    + I_neminus1_npminus1^2 * pz_plus * e_minus
                                    + 2 * I_ne_np * I_neminus1_npminus1 * (1+k_zp/(M_L_p_star+M_p_star))
                                        * (2*ELEM_CHARGE*B*sqrt(n_e*n_p))/((M_L_p_star+M_p_star)*E_Fe)))
    
        gv_minus_ga_term = (8 * (G_V - G_A)^2 * ((M_L_n_star - k_zn) * ((I_ne_np^2 * pB * e_plus)
                                    + (I_ne_npminus1^2 * pz_minus * e_plus))
                                    + (M_L_n_star + k_zn) * ((I_neminus1_np^2 * pB * e_minus)
                                    + (I_neminus1_npminus1^2 * pz_minus * e_minus))))

        gv2_minus_ga2_term = (16 * (G_V^2-G_A^2) * M_p_star * (I_ne_np^2 * pB * e_plus
                                    + I_neminus1_npminus1^2 * (k_zp^2/(M_L_p_star + M_p_star)^2 - 1) * e_minus
                                    + 2 * I_ne_np * I_neminus1_npminus1 * k_zp
                                        * (2*ELEM_CHARGE*B*sqrt(n_e*n_p))/((M_L_p_star+M_p_star)^2*E_Fe)))
    else
        @error "m_red: Proton spin not set to ±0.5."
    end

    gv_plus_ga_term + gv_minus_ga_term + gv2_minus_ga2_term
end


"""
    _sign_sum(spin, n_p, n_e, E_Fn, E_Fe, k_Fn, k_zp, k_ze, B, M_n_star, M_p_star)

Calculate the sum of the reduced matrix element over the signs of the proton and electron z-momenta.

Includes a catch to remove sign combinations that would violate conservation of momentum.
"""
function _sign_sum(spin, n_p, n_e, E_Fn, E_Fe, k_Fn, k_zp, k_ze, B, M_n_star, M_p_star)
    sum_over_signs = 0
    if k_Fn - abs(k_ze+k_zp) >= 0 # momentum-conserving step function -- abs() makes kze+kzp and -kzp-kze equivalent
        sum_over_signs += (m_red(spin, n_p, n_e, E_Fn, E_Fe, k_Fn, k_zp, k_ze, B, M_n_star, M_p_star)
            + m_red(spin, n_p, n_e, E_Fn, E_Fe, k_Fn, -k_zp, -k_ze, B, M_n_star, M_p_star))
    end
        
    if k_Fn - abs(k_ze-k_zp) >= 0 # momentum-conserving step function -- abs() makes kze-kzp and kzp-kze equivalent
        sum_over_signs += (m_red(spin, n_p, n_e, E_Fn, E_Fe, k_Fn, -k_zp, k_ze, B, M_n_star, M_p_star)
            + m_red(spin, n_p, n_e, E_Fn,  E_Fe, k_Fn, k_zp, -k_ze, B, M_n_star, M_p_star))
    end
    sum_over_signs
end


"""
    _integrand_func(k_zp_tilde, k_ze_tilde, n_p, n_e, k_Fp, k_Fe, M_p_star, proton_spin, B, T)

Calculate the integrand function for the phasespace integrals, ``-6*\\text{Li}_4(-\\exp(-x_e-x_p))*f(-x_e)*f(-x_p)``.
"""
function _integrand_func(k_zp_tilde, k_ze_tilde, n_p, n_e, k_Fp, k_Fe, M_p_star, proton_spin, B, T)
    x_e_val = sqrt(k_ze_tilde^2 + 2 * n_e * ELEM_CHARGE * B / T^2) - k_Fe / T
    x_p_val = sqrt(k_zp_tilde^2 + (M_p_star^2 + 2*n_p*ELEM_CHARGE*B - proton_spin*PROTON_G_MINUS_2*ELEM_CHARGE*B) / T^2) - (sqrt(k_Fp^2 + M_p_star^2) / T)

    if x_e_val < -150 || x_p_val < -150 || x_e_val + x_p_val > 150 # prevents OverflowError
        return 0
    end

    x_p_val > 150 ? fd_p = 1 : fd_p = 1 + exp(-x_p_val)
    x_e_val > 150 ? fd_e = 1 : fd_e = 1 + exp(-x_e_val)
    return -6 * reli4(-exp(-x_e_val - x_p_val)) / (fd_p * fd_e)
end


"""
    _integrand_func(k_zp_tilde, k_ze_tilde, n_p, n_e, k_Fp, k_Fe, M_p_star, proton_spin, B, T)

Calculate the integrand function for the phasespace integrals, including the sum over signs: ``-6*\\text{Li}_4(-\\exp(-x_e-x_p))*f(-x_e)*f(-x_p) \\Sigma``.
"""
function _integrand_func_with_m_red(k_zp_tilde, k_ze_tilde, n_p, n_e, k_Fn, k_Fp, k_Fe, E_Fn, M_n_star, M_p_star, proton_spin, B, T)
    x_e_val = sqrt(k_ze_tilde^2 + 2 * n_e * ELEM_CHARGE * B / T^2) - k_Fe / T
    x_p_val = sqrt(k_zp_tilde^2 + (M_p_star^2 + 2*n_p*ELEM_CHARGE*B - proton_spin*PROTON_G_MINUS_2*ELEM_CHARGE*B) / T^2) - (sqrt(k_Fp^2 + M_p_star^2) / T)

    if x_e_val < -150 || x_p_val < -150 || x_e_val + x_p_val > 150 # prevents OverflowError
        return 0
    end
    
    # to get rid of fermi_dirac() overhead, it's recreated here.
    # we don't need the other asymptotic condition because it'll be caught above
    x_p_val > 100 ? fd_p = 1 : fd_p = 1 + exp(-x_p_val)
    x_e_val > 100 ? fd_e = 1 : fd_e = 1 + exp(-x_e_val)
    return (-6 * reli4(-exp(-x_e_val - x_p_val)) / (fd_p * fd_e)
        * _sign_sum(proton_spin, n_p, n_e, E_Fn, k_Fe, k_Fn, k_zp_tilde*T, k_ze_tilde*T, B, M_n_star, M_p_star))
end


"""
    _phasespace_integrals(n_p, n_e, k_Fp, proton_spin, µ_e, M_p_star, B, T)

Compute the phasespace integrals, without M_red.
"""
function _phasespace_integrals(n_p, n_e, k_Fp, proton_spin, µ_e, M_p_star, B, T)
    kzp_tilde_bounds = bound_kz_tilde(k_Fp, M_p_star, n_p, proton_spin, B, T)
    kze_tilde_bounds = bound_kz_tilde(µ_e, 0, n_e, 0, B, T)
    
    integrand_func_wrapper(kz_tilde_vec) = _integrand_func(kz_tilde_vec[1], kz_tilde_vec[2], n_p, n_e, k_Fp, µ_e, M_p_star, proton_spin, B, T)

    integral_result, _ = hcubature(integrand_func_wrapper, (kzp_tilde_bounds[1], kze_tilde_bounds[1]), 
        (kzp_tilde_bounds[2], kze_tilde_bounds[2]), maxevals=2000)
    integral_result
end


"""
    _phasespace_integrals_with_m_red(n_p, n_e, k_Fn, k_Fp, E_Fn, proton_spin, µ_e, M_n_star, M_p_star, B, T)

Compute the phasespace integrals, with M_red inside the integral.
"""
function _phasespace_integrals_with_m_red(n_p, n_e, k_Fn, k_Fp, E_Fn, proton_spin, µ_e, M_n_star, M_p_star, B, T)
    kzp_tilde_bounds = bound_kz_tilde(k_Fp, M_p_star, n_p, proton_spin, B, T)
    kze_tilde_bounds = bound_kz_tilde(µ_e, 0, n_e, 0, B, T)

    integrand_func_wrapper(kz_tilde_vec) = _integrand_func_with_m_red(kz_tilde_vec[1], kz_tilde_vec[2],
                                                n_p, n_e, k_Fn, k_Fp, µ_e, E_Fn, M_n_star, M_p_star, proton_spin, B, T)

    integral_result, _ = hcubature(integrand_func_wrapper, (kzp_tilde_bounds[1], kze_tilde_bounds[1]), 
        (kzp_tilde_bounds[2], kze_tilde_bounds[2]); maxevals=500) # maxevals brought pretty low for speed
    integral_result
end

density(k_Fn, k_Fp) = (k_Fn^3 + k_Fp^3)/(3*pi^2*197.3^3)


"""
    _emissivity_singlepair_wrapper(spin, n_p, n_e, E_Fn, µ_e, k_Fp, k_Fn, B_MeV2, M_n_star, M_p_star, T, analytical_approx, m_red_in_int)

Compute the emissivity at a single LL pair ``(n_p, n_e)``.
"""
function _emissivity_singlepair_wrapper(spin, n_p, n_e, E_Fn, µ_e, k_Fp, k_Fn, B_MeV2, M_n_star, M_p_star, T, analytical_approx, m_red_in_int, force_numerical)
    k_zp = sqrt(max(0, k_Fp^2 - 2*n_p*ELEM_CHARGE*B_MeV2 + spin*PROTON_G_MINUS_2*ELEM_CHARGE*B_MeV2))
    k_ze = sqrt(max(0, µ_e^2 - 2*n_e*ELEM_CHARGE*B_MeV2))
    sum_over_signs = _sign_sum(spin, n_p, n_e, E_Fn, µ_e, k_Fn, k_zp, k_ze, B_MeV2, M_n_star, M_p_star)
    # Calculate Landau effective mass of the proton
    M_L_n_star = sqrt(M_n_star^2 + k_Fn^2)
    M_L_p_star = sqrt(M_p_star^2 + k_zp^2 + 2*n_p*ELEM_CHARGE*B_MeV2 - ELEM_CHARGE*B_MeV2*spin*PROTON_G_MINUS_2)
    # If the integral bounds extend to zero or k_z is small, need the numerical integration; otherwise, analytical is fine
    if force_numerical || !analytical_approx || k_ze/µ_e < 0.2 || k_zp/k_Fp < 0.2
        if m_red_in_int
            emissivity_from_current_pair = (NUMERICAL_PREFACTOR * (M_L_p_star+M_p_star)/M_L_p_star 
            * _phasespace_integrals_with_m_red(n_p, n_e, k_Fn, k_Fp, E_Fn, spin, µ_e, M_n_star, M_p_star, B_MeV2, T))
        else
            emissivity_from_current_pair = (NUMERICAL_PREFACTOR * (M_L_p_star+M_p_star)/M_L_p_star
            * _phasespace_integrals(n_p, n_e, k_Fp, spin, µ_e, M_p_star, B_MeV2, T) * sum_over_signs)
        end
    else
        emissivity_from_current_pair = ANALYTICAL_PREFACTOR * µ_e * (M_L_p_star + M_p_star) / abs(k_zp*k_ze) * sum_over_signs
    end
    emissivity_from_current_pair
end


"""
    direct_urca_rate(µ_B, k_Fn, k_Fp, µ_e, M_n_star, M_p_star, B, T;
                        spin_split=true, thermal_population=true, analytical_approx=true, m_red_in_int=true, cgs=true)

Calculate the emissivity of the Direct Urca process.
"""
function direct_urca_rate(µ_B, k_Fn, k_Fp, µ_e, M_n_star, M_p_star, B, T;
                        spin_split=true, thermal_population=true, analytical_approx=false, m_red_in_int=true, cgs=true) # take in B in -> GAUSS <-
    # Determine physical values based on the chosen toggles
    E_Fn = µ_B
    spin_split ? spin_up = 0.5 : spin_up = 0
    thermal_population ? therm_pop_temp = T : therm_pop_temp = 0

    B_MeV2 = B * GAUSS_TO_MEV2
    landau_sum = 0
    
    # Maximum Landau level attainable for proton (spin up/down) and electron
    # If spin splitting is disabled, n_max_p_up = n_max_p_down, and we can use
    # either one of them
    n_max_p_up = available_landau_levels(spin_up, k_Fp, M_p_star, B_MeV2, therm_pop_temp)
    n_max_p_down = available_landau_levels(-spin_up, k_Fp, M_p_star, B_MeV2, therm_pop_temp)
    n_max_e = available_landau_levels(0, µ_e, 0, B_MeV2, therm_pop_temp)
    
    if spin_split
        n_max_arr = [n_max_p_up n_max_p_down n_max_e]
    else
        n_max_arr = [n_max_p_up n_max_e]
    end
    
    # Sum up the contributions from each pair of LLs
    for n_e in 0:n_max_e
        for n_p_up in 0:n_max_p_up
            landau_sum += _emissivity_singlepair_wrapper(spin_up, n_p_up, n_e, E_Fn,
            µ_e, k_Fp, k_Fn, B_MeV2, M_n_star, M_p_star, T, analytical_approx, m_red_in_int, n_e >= n_max_e-1 || n_p_up >= n_max_p_up-1)
        end
        if spin_split
            for n_p_down in 0:n_max_p_down
                landau_sum += _emissivity_singlepair_wrapper(-spin_up, n_p_down, n_e, E_Fn,
                µ_e, k_Fp, k_Fn, B_MeV2, M_n_star, M_p_star, T, analytical_approx, m_red_in_int, n_e >= n_max_e-1 || n_p_down >= n_max_p_down-1)
            end
        end
    end
    
    if cgs
        MEV4_TO_CGS * (B_MeV2*T^6*landau_sum), n_max_arr
    else
        B_MeV2*T^6*landau_sum, n_max_arr
    end
end


"""
    find_x_from_µ(µ_B, B, mageos, nuc_inter)

Convert a baryon chemical potential ``µ_B`` to the corresponding value of Baiko and Yakovlev's parameter ``x``.
"""
function find_x_from_µ(µ_B, B, mageos, nuc_inter)
    if mageos
        eos_data = eos_data_mag(iufsu_star_constants_mag, µ_B/197.3, B/1e15) * 197.3
    else
        eos_data = eos_data_rmf(iufsu_star_constants, µ_B/197.3) * 197.3
    end

    # process EoS params 
    if nuc_inter
        k_Fn, k_Fp, µ_e, M_n_star, M_p_star = process_eos_data(eos_data, mageos)
        E_Fn = µ_B
        E_Fp = µ_B - µ_e 
    else
        k_Fn, k_Fp, µ_e, _, _ = process_eos_data(eos_data, mageos)
        M_n_star = 940.6
        M_p_star = 938.3
        E_Fn = µ_B
        E_Fp = µ_B - µ_e
    end

    BY_n = (k_Fp^2/(2*ELEM_CHARGE*B * GAUSS_TO_MEV2))
    (k_Fn^2 - (k_Fp+µ_e)^2)/(k_Fp^2 * BY_n^(-2/3))
end


"""
    find_µ_range_endpoints(x_min, x_max, B; mageos=true, nuc_inter=true)

Convert a baryon chemical potential ``µ_B`` to the corresponding value of Baiko and Yakovlev's parameter ``x``.
"""
function find_µ_range_endpoints(x_min, x_max, B; mageos=true, nuc_inter=true)
    µ1(µ) = find_x_from_µ(µ, B, mageos, nuc_inter) - x_min
    µ2(µ) = find_x_from_µ(µ, B, mageos, nuc_inter) - x_max
    (find_zero(µ2, 1100), find_zero(µ1, 1100))
end


"""
    rb(µ_B, B, T; toys=false, analytical_approx=false, m_red_in_int=false)

Find the parameters ``RB_\text{qc}`` and ``RB``. (Also returns ``x`` and ``\rho`` for convenience.)
"""
function rb(µ_B, B, T; toys=false, analytical_approx=false, m_red_in_int=false)
    # process EoS params 
    if toys
        eos_data = eos_data_mag(iufsu_star_constants_mag, µ_B/197.3, B/1e15) * 197.3
        k_Fn, k_Fp, µ_e, M_n_star, M_p_star = process_eos_data(eos_data, true)
        E_Fn = µ_B
        E_Fp = µ_B - µ_e 
    else
        eos_data = eos_data_rmf(iufsu_star_constants, µ_B/197.3) * 197.3
        k_Fn, k_Fp, µ_e, _, _ = process_eos_data(eos_data, false)
        M_n_star = NEUTRON_MASS
        M_p_star = PROTON_MASS
        E_Fn = µ_B
        E_Fp = µ_B - µ_e
    end

    BY_n = (k_Fp^2/(2*ELEM_CHARGE*B * GAUSS_TO_MEV2))
    BY_x = (k_Fn^2 - (k_Fp+µ_e)^2)/(k_Fp^2 * BY_n^(-2/3))
    mn_star = sqrt(k_Fn^2 + M_n_star^2) # *Landau* effective mass
    mp_star = sqrt(k_Fp^2 + M_p_star^2) # *Landau* effective mass

    Q_nu_0 = 457/10080 * pi * (G_F*cos(CABIBBO_ANGLE))^2 * (1+3*G_A^2) * mn_star * mp_star * µ_e * T^6

    # Our calculation
    if toys
        Q_nu, n_max_arr  = direct_urca_rate(µ_B, k_Fn, k_Fp, µ_e, M_n_star, M_p_star, B, T;
        spin_split=true, thermal_population=true, analytical_approx=analytical_approx, m_red_in_int=m_red_in_int, cgs=false)
    else
        Q_nu, n_max_arr = direct_urca_rate(µ_B, k_Fn, k_Fp, µ_e, M_n_star, M_p_star, B, T;
        spin_split=false, thermal_population=false, analytical_approx=analytical_approx, m_red_in_int=m_red_in_int, cgs=false)
    end

    # BY's calculation
    # Forbidden region
    if 0 <= BY_x
        airy_integrand(theta_s_vec) = sin(theta_s_vec[1])^(2/3) * airyai( (BY_x+theta_s_vec[2]^2) / (2^(4/3)*(sin(theta_s_vec[1]))^(2/3)) )^2
        rb_qc = 2^(-2/3) * hcubature(airy_integrand, (0, -100), (pi, 100))[1]
    # Permitted region
    else
        BY_phi = (1.211 + 0.4823*abs(BY_x) + 0.8453*abs(BY_x)^2.533) / (1 + 1.438*abs(BY_x)^1.209)
        rb_qc = 1 - cos(BY_phi)/((0.5816 + abs(BY_x))^1.192)
    end

    Q_nu/Q_nu_0, rb_qc, BY_x, density(k_Fn, k_Fp)
end


"""
    qc_emissivity(k_Fn, k_Fp, B, T; µ_e=-1)

Calculate the DUrca emissivity according to the quasiclassical formula from BY.
"""
function qc_emissivity(k_Fn, k_Fp, M_n_star, M_p_star, B, T; µ_e=-1)
    if µ_e == -1
        µ_e = k_Fp
    end
    BY_n = (k_Fp^2/(2*ELEM_CHARGE*B * GAUSS_TO_MEV2))
    BY_x = (k_Fn^2 - (k_Fp+µ_e)^2)/(k_Fp^2 * BY_n^(-2/3))
    mn_star = sqrt(k_Fn^2 + M_n_star^2) # *Landau* effective mass
    mp_star = sqrt(k_Fp^2 + M_p_star^2) # *Landau* effective mass
    
    Q_nu_0 = 457/10080 * pi * (G_F*cos(CABIBBO_ANGLE))^2 * (1+3*G_A^2) * mn_star * mp_star * µ_e * T^6
    
    # Forbidden region
    if 0 <= BY_x
        airy_integrand(theta_s_vec) = sin(theta_s_vec[1])^(2/3) * airyai( (BY_x+theta_s_vec[2]^2) / (2^(4/3)*(sin(theta_s_vec[1]))^(2/3)) )^2
        R_b = 2^(-2/3) * hcubature(airy_integrand, (0, -100), (pi, 100))[1]
    # Permitted region
    else
        BY_phi = (1.211 + 0.4823*abs(BY_x) + 0.8453*abs(BY_x)^2.533) / (1 + 1.438*abs(BY_x)^1.209)
        R_b = 1 - cos(BY_phi)/((0.5816 + abs(BY_x))^1.192)
    end
    
    (Q_nu_0 * R_b * MEV4_TO_CGS), R_b, BY_x
end


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


"""
    urca_rate_wrapper(mu_B, B, T; mageos=true, nuc_inter=true, spin_split=true, thermal_population=true, analytical_approx=false, m_red_in_int=true)

Compute the direct Urca, modified Urca, and quasiclassical DUrca rates, as well as the density and number of available LLs, for a given baryon chemical potential.
"""
function urca_rate_wrapper(mu_B, B, T; mageos=true, nuc_inter=true, spin_split=true, thermal_population=true, analytical_approx=false, m_red_in_int=true)
    # calculate EoS params
    if mageos
        eos_data = eos_data_mag(iufsu_star_constants_mag, mu_B/197.3, B/1e15) * 197.3
    else
        eos_data = eos_data_rmf(iufsu_star_constants, mu_B/197.3) * 197.3
    end

    # process EoS params 
    if nuc_inter
        k_Fn, k_Fp, µ_e, M_n_star, M_p_star = process_eos_data(eos_data, mageos)
        E_Fn = mu_B
        E_Fp = mu_B - µ_e 
    else
        k_Fn, k_Fp, µ_e, _, _ = process_eos_data(eos_data, mageos)
        M_n_star = 940.6
        M_p_star = 938.3
        E_Fn = mu_B
        E_Fp = mu_B - µ_e
    end
    density_val = density(k_Fn, k_Fp)

    # calculate rate values for both reactions
    DUrca_rate_val, n_maxes = direct_urca_rate(mu_B, k_Fn, k_Fp, µ_e, M_n_star, M_p_star, B, T;
        spin_split=spin_split, thermal_population=thermal_population, analytical_approx=analytical_approx, m_red_in_int=m_red_in_int, cgs=true)
    MUrca_rate_val = modified_urca_rate(k_Fn, k_Fp, M_n_star, M_p_star, T, cgs=true)
    qc_rate_val = qc_emissivity(k_Fn, k_Fp, M_n_star, M_p_star, B, T, µ_e=µ_e)[1]

    # println("durca rate args: ", mu_B, " ", k_Fn, " ", k_Fp, " ", µ_e, " ", M_n_star, " ", M_p_star)
    DUrca_rate_val, MUrca_rate_val, qc_rate_val, density_val, n_maxes #, durca_threshold_val
end
end # module Urca