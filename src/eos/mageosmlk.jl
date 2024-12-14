using NLsolve: nlsolve

#NEXT SECTION IS CONSTANTS AND BASIC FUNCTIONS

#all quantities in fm unless otherwise noted
NEUTRON_MASS = 939.57/197.3
PROTON_MASS = 938.27/197.3
NUCLEON_MASS = (NEUTRON_MASS + PROTON_MASS) / 2
MUON_MASS = 105.7/197.3
ELEM_CHARGE = sqrt(4 * pi / 137)
NEUTRON_MAGNETIC_MOMENT = -3.826 * ELEM_CHARGE / (2 * NEUTRON_MASS)
PROTON_G_MINUS_2 = 3.586
#conversion from 10^15 G to fm
G15_TOFM2 = 5.01e-4

#maximum available landau level accounting for g-2
function _max_landau_level(bare_energy, mass, b, g_minus_2, spin)
    if bare_energy^2 - mass^2 + spin * ELEM_CHARGE * b * g_minus_2 / 2 < 0
        0
    else
        floor(Int, (bare_energy^2 - mass^2 + spin * ELEM_CHARGE * b * g_minus_2 / 2)/(2 * ELEM_CHARGE *b))
    end
end

#density of a neutral spin half fermion accounting for magnetic moment
function _neutral_density(bare_energy, mass, b, g_mu)
    if g_mu > 0
        if bare_energy - g_mu * b / 2 < mass
            ((bare_energy + g_mu * b / 2)^2 - mass^2)^(3/2) / (6*pi^2)
        else
            ((bare_energy + g_mu * b / 2)^2 - mass^2)^(3/2) / (6*pi^2) + ((bare_energy - g_mu * b / 2)^2 - mass^2)^(3/2) / (6*pi^2)
        end
    else
        if bare_energy + g_mu * b / 2 < mass
            ((bare_energy + g_mu * b / 2)^2 - mass^2)^(3/2) / (6*pi^2)
        else
            ((bare_energy + g_mu * b / 2)^2 - mass^2)^(3/2) / (6*pi^2) + ((bare_energy - g_mu * b / 2)^2 - mass^2)^(3/2) / (6*pi^2)
        end
    end
end

#density of a charged spin half fermion, including g-2 and landau levels
function _charged_density(bare_energy, mass, b, g_minus_2)
    #method to include g-2 only works for protons, set g_minus_2 to zero for electrons and muons
    val = 0
    anom_mom_energy = g_minus_2 * ELEM_CHARGE * b / 2
    
    if bare_energy^2 + anom_mom_energy > mass^2
        #zeroth Landau level does not have spin degeneracy
        kz_temp = sqrt(bare_energy^2 - mass^2 + anom_mom_energy)
        #extra factor of 2 because kz can be positive or negative
        val += ELEM_CHARGE * b / (4*pi^2) * 2 * kz_temp
    end

    for spin in (-1, 1)
        #highest available landau level
        n_max = _max_landau_level(bare_energy, mass, b, g_minus_2, spin)
        for n in range(1, n_max)
            kz_temp = sqrt(bare_energy^2 - mass^2 - 2 * n * ELEM_CHARGE * b + spin * anom_mom_energy)
            val += ELEM_CHARGE * b / (4 * pi^2) * 2 * kz_temp
        end
    end
    val
end

#scalar density of a neutral spin half fermion, accounting for magnetic moment
function _neutral_scalar_density(bare_energy, mass, b, g_mu)
    #different fermi momenta for different spin orientation
    if g_mu > 0
        kf_up = sqrt((bare_energy + g_mu * b/2)^2 - mass^2)
        if kf_up / mass > 0.1
            scalar_density_up = mass / (4*pi^2) * (kf_up * sqrt(kf_up^2 + mass^2) - mass^2 * atanh(kf_up / sqrt(kf_up^2 + mass^2)))
        else
            scalar_density_up = kf_up^3 / (3 * pi^2) - kf_up^5 / (10 * pi^2 * mass^2)
        end

        if bare_energy - g_mu * b / 2 < mass
            scalar_density_up
        else
            kf_down = sqrt((bare_energy - g_mu * b/2)^2 - mass^2)
            if kf_down / mass > 0.1
                scalar_density_down = mass / (4*pi^2) * (kf_down * sqrt(kf_down^2+mass^2) - mass^2 * atanh(kf_down / sqrt(kf_down^2 + mass^2)))
            else
                scalar_density_down = kf_down^3 / (3 * pi^2) - kf_down^5 / (10 * pi^2 * mass^2)
            end
            scalar_density_up + scalar_density_down
        end
    else
        kf_down = sqrt((bare_energy - g_mu * b/2)^2 - mass^2)
        if kf_down / mass > 0.1
            scalar_density_down = mass / (4*pi^2) * (kf_down * sqrt(kf_down^2+mass^2) - mass^2 * atanh(kf_down / sqrt(kf_down^2 + mass^2)))
        else
            scalar_density_down = kf_down^3 / (3 * pi^2) - kf_down^5 / (10 * pi^2 * mass^2)
        end

        if bare_energy + g_mu * b / 2 < mass
            scalar_density_up
        else
            kf_up = sqrt((bare_energy + g_mu * b/2)^2 - mass^2)
            if kf_up / mass > 0.1
                scalar_density_up = mass / (4*pi^2) * (kf_up * sqrt(kf_up^2 + mass^2) - mass^2 * atanh(kf_up / sqrt(kf_up^2 + mass^2)))
            else
                scalar_density_up = kf_up^3 / (3 * pi^2) - kf_up^5 / (10 * pi^2 * mass^2)
            end
            scalar_density_up + scalar_density_down
        end
    end
end

#scalar density fo a charged spin half fermion, account for landau levels and g-2
function _charged_scalar_density(bare_energy, mass, b, g_minus_2)
    val=0
    anom_mom_energy = g_minus_2 * ELEM_CHARGE * b / 2

    if bare_energy^2 + g_minus_2 * ELEM_CHARGE * b / 2 > mass^2
        #zeroth Landau level does not have spin degeneracy
        kz_temp = sqrt(bare_energy^2 - mass^2 + anom_mom_energy)
        if kz_temp / sqrt(mass^2 - anom_mom_energy) > 0.1
            #extra factor of 2 since kz can be positive or negative
            val += ELEM_CHARGE * b / (4*pi^2) * mass * 2 * asinh(kz_temp / sqrt(mass^2 - anom_mom_energy))
        else
            val += ELEM_CHARGE * b / (4 * pi^2) * mass * 2 * kz_temp / sqrt(mass^2 - anom_mom_energy)
        end
    end

    for spin in (-1, 1)
        #highest available landau level
        n_max = _max_landau_level(bare_energy, mass, b, g_minus_2, spin)
        for n in range(1, n_max)
            kz_temp = sqrt(bare_energy^2 - mass^2 - 2 * n * ELEM_CHARGE * b + spin * anom_mom_energy)
            if kz_temp / sqrt(mass^2 + 2 * n * ELEM_CHARGE * b - spin * g_minus_2 * ELEM_CHARGE * b) < 0.1
                val += ELEM_CHARGE * b / (4 * pi^2) * mass * 2 * kz_temp / sqrt(mass^2 
                    + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy)
            else
                val += ELEM_CHARGE * b / (4 * pi^2) * mass * 2 * asinh(kz_temp / sqrt(mass^2 + 2 * n * ELEM_CHARGE * b 
                    - spin * anom_mom_energy))
            end
        end
    end
    val
end

#energy density of neutral spin half fermion accounting for magnetic moment
function _neutral_energy_density(bare_energy, mass, b, g_mu)
    if g_mu > 0
        kf_up = sqrt((bare_energy + g_mu * b/2)^2 - mass^2)
        if kf_up / mass > 0.1
            energy_density_up = 1/(8*pi^2) * (kf_up * sqrt(kf_up^2 + mass^2) * (2*kf_up^2 + mass^2) - mass^4 * asinh(kf_up/mass))
        else
            energy_density_up = mass * kf_up^3 / (6 * pi^2) + 1 / (10 * pi^2) * kf_up^5 / mass
        end

        if bare_energy - g_mu * b / 2 < mass
            energy_density_up
        else
            kf_down = sqrt((bare_energy - g_mu * b/2)^2 - mass^2)
            if kf_down / mass > 0.1
                energy_density_down = 1/(8*pi^2) * (kf_down * sqrt(kf_down^2 + mass^2) * (2*kf_down^2 + mass^2) - mass^4 * asinh(kf_down/mass))
            else
                energy_density_down = mass * kf_down^3 / (6 * pi^2) + 1 / (10 * pi^2) * kf_down^5 / mass
            end
            energy_density_up + energy_density_down
        end
    else
        kf_down = sqrt((bare_energy - g_mu * b/2)^2 - mass^2)
        if kf_down / mass > 0.1
            energy_density_down = 1/(8*pi^2) * (kf_down * sqrt(kf_down^2 + mass^2) * (2*kf_down^2 + mass^2) - mass^4 * asinh(kf_down/mass))
        else
            energy_density_down = mass * kf_down^3 / (6 * pi^2) + 1 / (10 * pi^2) * kf_down^5 / mass
        end

        if bare_energy + g_mu * b / 2 < mass
            energy_density_down
        else
            kf_up = sqrt((bare_energy + g_mu * b/2)^2 - mass^2)
            if kf_up / mass > 0.1
                energy_density_up = 1/(8*pi^2) * (kf_up * sqrt(kf_up^2 + mass^2) * (2*kf_up^2 + mass^2) - mass^4 * asinh(kf_up/mass))
            else
                energy_density_up = mass * kf_up^3 / (6 * pi^2) + 1 / (10 * pi^2) * kf_up^5 / mass
            end
            energy_density_up + energy_density_down
        end
    end
end

#energy density of a charged spin half fermion accounting for landau levels and g-2
function _charged_energy_density(bare_energy, mass, b, g_minus_2)
    val = 0
    anom_mom_energy = g_minus_2 * ELEM_CHARGE * b /2
    if bare_energy^2 + anom_mom_energy > mass^2
        #zeroth Landau level does not have spin degeneracy
        kz_temp = sqrt(bare_energy^2 - mass^2 + anom_mom_energy)
        if kz_temp / sqrt(mass^2 - anom_mom_energy) > 0.1
            val += ELEM_CHARGE * b / (4*pi^2) * (sqrt(kz_temp^2 + mass^2 - anom_mom_energy) * kz_temp + 
                (mass^2 - anom_mom_energy) * asinh(kz_temp / sqrt(mass^2 - anom_mom_energy)))
        else
            val += ELEM_CHARGE * b / (4 * pi^2) * (sqrt(kz_temp^2 + mass^2 - anom_mom_energy) * kz_temp + 
                sqrt(mass^2 - anom_mom_energy) * kz_temp)
        end
    end
    
    for spin in (-1, 1)
        n_max = _max_landau_level(bare_energy, mass, b, g_minus_2, spin)
        for n in range(1, n_max)
            kz_temp = sqrt(bare_energy^2 - mass^2 - 2 * n * ELEM_CHARGE * b + spin * anom_mom_energy)
            if kz_temp / sqrt(mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) > 0.1
                val+= ELEM_CHARGE * b / (4 * pi^2) * (sqrt(kz_temp^2 + mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * kz_temp 
                    + (mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * 
                    asinh(kz_temp / sqrt(mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy)))
            else
                val += ELEM_CHARGE * b / (4 * pi^2) * (sqrt(kz_temp^2 + mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * kz_temp
                    + sqrt(mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * kz_temp)
            end
        end
    end
    
    val
end

#pressure due to a neutral spin half fermion accounting for magnetic moment
function _neutral_pressure(bare_energy, mass, b, g_mu)
    if g_mu > 0
        kf_up = sqrt((bare_energy + g_mu * b/2)^2 - mass^2)
        if kf_up / mass > 0.1
            pressure_up = 1/(48*pi^2) * (kf_up * (2 * kf_up^2 - 3 * mass^2) * sqrt(kf_up^2 + mass^2) + (
                3 * mass^4 * atanh(kf_up / sqrt(kf_up^2 + mass^2))))
        else
            pressure_up = 1 / (15 * pi^2) * kf_up^5 / mass
        end

        if bare_energy - g_mu * b / 2 < mass
            pressure_up
        else
            kf_down = sqrt((bare_energy - g_mu * b/2)^2 - mass^2)
            if kf_down / mass > 0.1
                pressure_down = 1/(48*pi^2) * (kf_down * (2 * kf_down^2 - 3 * mass^2) * sqrt(kf_down^2 + mass^2) + (
                    3 * mass^4 * atanh(kf_down / sqrt(kf_down^2 + mass^2))))
            else
                pressure_down = 1 / (15 * pi^2) * kf_down^5 / mass
            end
            pressure_up + pressure_down
        end
    else
        kf_down = sqrt((bare_energy - g_mu * b/2)^2 - mass^2)
        if kf_down / mass > 0.1
            pressure_down = 1/(48*pi^2) * (kf_down * (2 * kf_down^2 - 3 * mass^2) * sqrt(kf_down^2 + mass^2) + (
                3 * mass^4 * atanh(kf_down / sqrt(kf_down^2 + mass^2))))
        else
            pressure_down = 1 / (15 * pi^2) * kf_down^5 / mass 
        end
    
        if bare_energy + g_mu * b / 2 < mass
            pressure_down
        else
            kf_up = sqrt((bare_energy + g_mu * b/2)^2 - mass^2)
            if kf_up / mass > 0.1
                pressure_up = 1/(48*pi^2) * (kf_up * (2 * kf_up^2 - 3 * mass^2) * sqrt(kf_up^2 + mass^2) + (
                    3 * mass^4 * atanh(kf_up / sqrt(kf_up^2 + mass^2))))
            else
                pressure_up = 1 / (15 * pi^2) * kf_up^5 / mass
            end
            pressure_up + pressure_down
        end
    end
end
    
#pressure from a spin half charged fermion account for landau levels and g-2
function _charged_pressure(bare_energy, mass, b, g_minus_2)
    val = 0

    anom_mom_energy = g_minus_2 * ELEM_CHARGE * b /2
    if bare_energy^2 + anom_mom_energy > mass^2
        #zeroth Landau level does not have spin degeneracy
        kz_temp = sqrt(bare_energy^2 - mass^2 + anom_mom_energy)
        if kz_temp / sqrt(mass^2 - anom_mom_energy) > 0.1
            val += ELEM_CHARGE * b / (12*pi^2) * (sqrt(kz_temp^2 + mass^2 - anom_mom_energy) * kz_temp
                - (mass^2 - anom_mom_energy) * asinh(kz_temp / sqrt(mass^2  - anom_mom_energy)))
        else
            val += ELEM_CHARGE * b / (12 * pi^2) * (sqrt(kz_temp^2 + mass^2 - anom_mom_energy) * kz_temp + 
                - sqrt(mass^2 - anom_mom_energy) * kz_temp)
        end
    end
    
    for spin in (-1, 1)
        n_max = _max_landau_level(bare_energy, mass, b, g_minus_2, spin)
        for n in range(1, n_max)
            kz_temp = sqrt(bare_energy^2 - mass^2 - 2 * n * ELEM_CHARGE * b + spin * anom_mom_energy)
            if kz_temp / sqrt(mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) > 0.1
                val+= ELEM_CHARGE * b / (12*pi^2) * (sqrt(kz_temp^2 + mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * kz_temp
                    - (mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * asinh(kz_temp / sqrt(mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy)))
            else
                val += ELEM_CHARGE * b / (4 * pi^2) * (sqrt(kz_temp^2 + mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * kz_temp
                    - sqrt(mass^2 + 2 * n * ELEM_CHARGE * b - spin * anom_mom_energy) * kz_temp)
            end
        end
    end
    val
end

#number density of leptons in a magnetic field with a given chemical potential
function _lepton_number_density(mue, b)
    val = 0
    ne_max = _max_landau_level(mue, 0, b, 0, 0)
    #zeroth Landau level does not have spin degeneracy
    kz_temp = mue
    #extra factor of two for signs of kz
    val += ELEM_CHARGE * b / (4*pi^2) * 2 * kz_temp
    #can add all the way up to n_max since electron g-2 is approx 0
    for n in range(1, ne_max)
        kz_temp = sqrt(mue^2 - 2 * n * ELEM_CHARGE * b)
        val += ELEM_CHARGE * b / (2*pi^2) * 2 * kz_temp
    end
    #if mue is larger than the muon mass, include muons as well
    if mue > MUON_MASS
        val += _charged_density(mue, MUON_MASS, b, 0)
    end
    val
end

#energy density of leptons in a magnetic field with a given chemical potential
function _lepton_energy_density(mue, b)
    val = 0
    ne_max = _max_landau_level(mue, 0, b, 0, 0)
    #zeroth Landau level does not have spin degeneracy
    kz_temp = mue
    val += ELEM_CHARGE * b / (4*pi^2) * kz_temp^2
    #can add all the way up to n_max since electron g-2 is approx 0
    for n in range(1, ne_max)
        kz_temp = sqrt(mue^2 - 2 * n * ELEM_CHARGE * b)
        if kz_temp / sqrt(2 * n * ELEM_CHARGE * b) > 0.1
            val += ELEM_CHARGE * b / (2*pi^2) * (sqrt(kz_temp^2 + 2 * n * ELEM_CHARGE * b) * kz_temp
                + (2 * n * ELEM_CHARGE * b) * asinh(kz_temp / sqrt(2 * n * ELEM_CHARGE * b)))
        else
            val += ELEM_CHARGE * b / (2 * pi^2) * (sqrt(kz_temp^2 + 2 * n * ELEM_CHARGE * b) * kz_temp
                + sqrt(2 * n * ELEM_CHARGE * b) * kz_temp)
        end
    end
    #if mue is greater than the muon mass, include muons
    if mue > MUON_MASS
        val += _charged_energy_density(mue, MUON_MASS, b, 0)
    end
    val
end

#pressure of leptons in a magnetic field with a given chemical potential
function _lepton_pressure(mue, b)
    val = 0
    ne_max = _max_landau_level(mue, 0, b, 0, 0)
    #zeroth Landau level does not have spin degeneracy
    kz_temp = mue
    val += ELEM_CHARGE * b / (12*pi^2) * kz_temp^2
    #can add all the way up to n_max since electron g-2 is approx 0
    for n in range(1, ne_max)
        kz_temp = sqrt(mue^2 - 2 * n * ELEM_CHARGE * b)
        if kz_temp / sqrt(2 * n * ELEM_CHARGE * b) > 0.1
            val += ELEM_CHARGE * b / (6*pi^2) * (sqrt(kz_temp^2 + 2 * n * ELEM_CHARGE * b) * kz_temp
                - (2 * n * ELEM_CHARGE * b) * asinh(kz_temp / sqrt(2 * n * ELEM_CHARGE * b)))
        else
            val += ELEM_CHARGE * b / (6 * pi^2) * (sqrt(kz_temp^2 + 2 * n * ELEM_CHARGE * b) * kz_temp 
                - sqrt(2 * n * ELEM_CHARGE * b) * kz_temp)
        end
    end
    if mue > MUON_MASS
        val += _charged_pressure(mue, MUON_MASS, b, 0)
    end
    val
end

function _rmf_energy_density_mag(fields, constants, b)
    #rmf constants
    msig,mw,mrho,gsig,gw,grho,kappa3,kappa4,zeta0,eta1,eta2,etarho,eta1rho,eta2rho = constants
    #meson and nucleon fields
    sigma, omega, rho, proton_bare_energy, neutron_bare_energy, mue = fields

    #nucleon effective masses
    mp_star = PROTON_MASS - gsig * sigma
    mn_star = NEUTRON_MASS - gsig * sigma

    #nucleon and meson energy densities
    nucleon_energy_density = _neutral_energy_density(neutron_bare_energy, mn_star, b, NEUTRON_MAGNETIC_MOMENT) + (
        _charged_energy_density(proton_bare_energy, mp_star, b, PROTON_G_MINUS_2))
    meson_energy_density = (1/2 * (msig^2*sigma^2 + mw^2*omega^2 + mrho^2*rho^2) + kappa3/(6*NUCLEON_MASS)*gsig*msig^2*sigma^3 + kappa4/(24*NUCLEON_MASS^2)*gsig^2*msig^2*sigma^4
        + zeta0/8*gw^2*omega^4 + eta1/(2*NUCLEON_MASS)*gsig*mw^2*sigma*omega^2 + eta2/(4*NUCLEON_MASS^2)*gsig^2*mw^2*sigma^2*omega^2
        + etarho/(2*NUCLEON_MASS)*gsig*mrho^2*sigma*rho^2 + eta1rho/(4*NUCLEON_MASS^2)*gsig^2*mrho^2*sigma^2*rho^2 + 3*eta2rho/(4*NUCLEON_MASS^2)*gw^2*mrho^2*omega^2*rho^2)

    #add energy density due to leptons
    nucleon_energy_density + meson_energy_density + _lepton_energy_density(mue, b)
end

function _rmf_pressure_mag(fields, constants, b)
    #rmf constants
    msig,mw,mrho,gsig,gw,grho,kappa3,kappa4,zeta0,eta1,eta2,etarho,eta1rho,eta2rho = constants
    #meson and nucleon fields
    sigma, omega, rho, proton_bare_energy, neutron_bare_energy, mue = fields

    #nucleon effective masses
    mp_star = PROTON_MASS - gsig * sigma
    mn_star = NEUTRON_MASS - gsig * sigma

    #nucleon and meson pressures
    nucleon_pressure = _neutral_pressure(neutron_bare_energy, mn_star, b, NEUTRON_MAGNETIC_MOMENT) + (
        _charged_pressure(proton_bare_energy, mp_star, b, PROTON_G_MINUS_2))
    meson_pressure = (1/2 * (-msig^2*sigma^2 + mw^2*omega^2 + mrho^2*rho^2) - kappa3/(6*NUCLEON_MASS)*gsig*msig^2*sigma^3 - kappa4/(24*NUCLEON_MASS^2)*gsig^2*msig^2*sigma^4
        + zeta0/24*gw^2*omega^4 + eta1/(2*NUCLEON_MASS)*gsig*mw^2*sigma*omega^2 + eta2/(4*NUCLEON_MASS^2)*gsig^2*mw^2*sigma^2*omega^2
        + etarho/(2*NUCLEON_MASS)*gsig*mrho^2*sigma*rho^2 + eta1rho/(4*NUCLEON_MASS^2)*gsig^2*mrho^2*sigma^2*rho^2 + eta2rho/(4*NUCLEON_MASS^2)*gw^2*mrho^2*omega^2*rho^2)
    
    #add pressure from leptons
    nucleon_pressure + meson_pressure + _lepton_pressure(mue, b)
end

function _rmf_eom_magnetic(fields, mub, constants, b)
    #rmf constants
    msig,mw,mrho,gsig,gw,grho,kappa3,kappa4,zeta0,eta1,eta2,etarho,eta1rho,eta2rho = constants
    #meson and nucleon fields
    sigma, omega, rho, proton_bare_energy, neutron_bare_energy, mue = fields
    
    #nucleon effective masses
    mp_star = PROTON_MASS - gsig * sigma
    mn_star = NEUTRON_MASS - gsig * sigma

    #throw out unphysical values of parameters
    if proton_bare_energy < mp_star || neutron_bare_energy < mn_star || mue < 0 || mp_star < 0 || mn_star < 0
        (10, 10, 10, 10, 10, 10)
    else
    #meson equations of motion
        val = [(gsig * (_charged_scalar_density(proton_bare_energy, mp_star, b, PROTON_G_MINUS_2) + _neutral_scalar_density(neutron_bare_energy, mn_star, b, NEUTRON_MAGNETIC_MOMENT))
            - msig^2*sigma - kappa3/(2*NUCLEON_MASS)*gsig*msig^2*sigma^2 - kappa4/(6*NUCLEON_MASS^2)*gsig^2*msig^2*sigma^3 
            + eta1/(2*NUCLEON_MASS)*gsig*mw^2*omega^2 + eta2/(2*NUCLEON_MASS^2)*gsig^2*mw^2*sigma*omega^2
            + etarho/(2*NUCLEON_MASS)*gsig*mrho^2*rho^2 + eta1rho/(2*NUCLEON_MASS^2)*gsig^2*mrho^2*sigma*rho^2)]
        push!(val, -gw * (_charged_density(proton_bare_energy, mp_star, b, PROTON_G_MINUS_2) + _neutral_density(neutron_bare_energy, mn_star, b, NEUTRON_MAGNETIC_MOMENT))
            + mw^2*omega + eta1/NUCLEON_MASS*gsig*mw^2*sigma*omega + eta2/(2*NUCLEON_MASS^2)*gsig^2*mw^2*sigma^2*omega
            + eta2rho/(2*NUCLEON_MASS^2)*gw^2*mrho^2*omega*rho^2 + 1/6*zeta0*gw^2*omega^3)
        push!(val, -1/2 * grho * (_charged_density(proton_bare_energy, mp_star, b, PROTON_G_MINUS_2) - _neutral_density(neutron_bare_energy, mn_star, b, NEUTRON_MAGNETIC_MOMENT))
            + mrho^2*rho + etarho/NUCLEON_MASS*gsig*mrho^2*sigma*rho + eta1rho/(2*NUCLEON_MASS^2)*gsig^2*mrho^2*sigma^2*rho
            + eta2rho/(2*NUCLEON_MASS^2)*gw^2*mrho^2*omega^2*rho)
        #add beta equilibrium
        push!(val, neutron_bare_energy + gw * omega - 1/2 * grho * rho - mub)
        push!(val, proton_bare_energy + gw * omega + 1/2 * grho * rho - mub + mue)
        #impose charge neutrality
        push!(val, _lepton_number_density(mue, b) - _charged_density(proton_bare_energy, mp_star, b, PROTON_G_MINUS_2))
        val
    end
end

struct RMFModelMagnetic
    #initialize RMF constants, does not include delta meson
    msig::Float64
    mw::Float64
    mrho::Float64
    gsig::Float64
    gw::Float64
    grho::Float64
    kappa3::Float64
    kappa4::Float64
    zeta0::Float64
    eta1::Float64
    eta2::Float64
    etarho::Float64
    eta1rho::Float64
    eta2rho::Float64
    ref::String
end
    
#solves for equation of state data for a given magnetic field and mub
function eos_data_mag(model::RMFModelMagnetic, mub, bg15, nuc_guess = [0.1, 0.2, -0.05, 5, 6, 1])
    #mag field should be in units of 10^15G, convert to fm
    b = G15_TOFM2 * bg15
    #list of constants to be passed to RMF functions
    constants = [model.msig,model.mw,model.mrho,model.gsig,model.gw,model.grho,model.kappa3,model.kappa4,model.zeta0,
                model.eta1,model.eta2,model.etarho,model.eta1rho,model.eta2rho]
    #solve for nucleon and meson fields
    fields_temp = nlsolve(x -> _rmf_eom_magnetic(x, mub, constants, b), nuc_guess)

    mp_star = PROTON_MASS - model.gsig * fields_temp.zero[1]
    mn_star = NEUTRON_MASS - model. gsig * fields_temp.zero[1]
    mue = fields_temp.zero[6] 
    proton_bare_energy = fields_temp.zero[4]
    neutron_bare_energy = fields_temp.zero[5]
    pressure = _rmf_pressure_mag(fields_temp.zero, constants, b)
    energy_density = _rmf_energy_density_mag(fields_temp.zero, constants, b)
    proton_density = _charged_density(proton_bare_energy, mp_star, b, PROTON_G_MINUS_2)
    neutron_density = _neutral_density(neutron_bare_energy, mn_star, b, NEUTRON_MAGNETIC_MOMENT)

    [proton_bare_energy, neutron_bare_energy, proton_density, neutron_density, mue, energy_density, pressure, mp_star, mn_star]
end

bsp_constants_mag = RMFModelMagnetic(NUCLEON_MASS*0.5383, NUCLEON_MASS*0.8333, NUCLEON_MASS*0.82, 4*pi*0.8764, 4*pi*1.1481,
    4*pi*1.0508, 1.0681, 14.9857, 0, 0.0872, 3.1265, 0, 0, 53.7642, "arXiv:1204.2644")

iufsu_star_constants_mag = RMFModelMagnetic(NUCLEON_MASS*0.543, NUCLEON_MASS*0.8331, NUCLEON_MASS*0.8198, 4*pi*0.8379, 4*pi*1.0666,
    4*pi*0.9889, 1.1418, 1.0328, 5.3895, 0, 0, 0, 0, 41.3066, "arXiv:1204.2644")

nl3_constants_mag = RMFModelMagnetic(NUCLEON_MASS*0.5412, NUCLEON_MASS*0.8333, NUCLEON_MASS*0.8126, 4*pi*0.8131, 4*pi*1.024,
    4*pi*0.7121, 1.4661, -5.6718, 0, 0, 0, 0, 0, 0, "arXiv:1204.2644")
        
