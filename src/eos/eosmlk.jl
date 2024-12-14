using NLsolve: nlsolve

#all quantities are in units of fm.  Use 1=hbar*c=197.3 MeV*fm to convert as necessary when finished

#NEXT SECTION IS FUNCTIONS AND CONSTANTS NEEDED FOR ALL MODELS

NEUTRON_MASS = 939.57 / 197.3
PROTON_MASS = 938.27 / 197.3
NUCLEON_MASS = (NEUTRON_MASS + PROTON_MASS) / 2
MUON_MASS = 105.7 / 197.3

#calculates the number density of electrons and muons

_density(kf) = kf^3 / (3 * pi^2)

function _lepton_density(mue)
    if mue > MUON_MASS
        _density(mue) + _density(sqrt(mue^2 - MUON_MASS^2))
    else
        _density(mue)
    end
end

function _fermion_energy_density(kf, mass) #includes spin 1/2 degeneracy but not color/flavor
    if kf / mass > 0.1
        1 / (8 * pi^2) * (kf * sqrt(kf^2 + mass^2) * (2 * kf^2 + mass^2) - mass^4 * asinh(kf / mass))
    else
        mass * _density(kf) + 1 / (10 * pi^2) * kf^5 / mass
    end 
end

function _lepton_energy_density(mue) #include muons if mue is larger than the muon mass
    if mue > MUON_MASS
        mue^4 / (4*pi^2) + _fermion_energy_density(sqrt(mue^2-MUON_MASS^2), MUON_MASS)
    else
        mue^4 / (4*pi^2)
    end
end

function _fermion_pressure(kf, mass) #includes spin 1/2 degeneracy, but not color/flavor degeneracy
    if kf / mass > 0.1
        1 / (24 * pi^2) * (kf * (2 * kf^2 - 3 * mass^2) * sqrt(kf^2 + mass^2) + 3 * mass^4 * asinh(kf / mass))
    else
        1 / (15 * pi^2) * kf^5 / mass
    end
end

function _lepton_pressure(mue) #include muons if mue is larger than the muon mass
    if mue > MUON_MASS
        mue^4 / (12 * pi^2) + _fermion_pressure(sqrt(mue^2-MUON_MASS^2), MUON_MASS)
    else
        mue^4 / (12 * pi^2)
    end
end
    
#scalar density needed for RMF equations of motion
function _scalar_density(kf, mass)
    if kf / mass > 0.1
        mass / (2 * pi^2) * (kf * sqrt(kf^2 + mass^2) - mass^2 * asinh(kf / mass))
    else
        _density(kf) - kf^5 / (10 * pi^2 * mass^2)
    end
end

function _rmf_eom_neutral(fields, mub, constants) #rmf eoms for charge neutral nuclear matter
    #rmf model constants
    msig, mw, mrho, gsig, gw, grho, kappa3, kappa4, zeta0, eta1, eta2, etarho, eta1rho, eta2rho = constants
    #fields to solve for
    sigma, omega, rho, kfp, kfn, mue = fields

    #effective nucleon masses
    mp_star = PROTON_MASS - sigma * gsig
    mn_star = NEUTRON_MASS - sigma * gsig

    #meson equations of motion
    val = [(gsig * (_scalar_density(kfp, mp_star) + _scalar_density(kfn, mn_star)) - msig^2 * sigma - kappa3 / (2 * NUCLEON_MASS) * gsig * msig^2 * sigma^2 
        - kappa4 / (6 * NUCLEON_MASS^2) * gsig^2 * msig^2 * sigma^3 + eta1 / (2 * NUCLEON_MASS) * gsig * mw^2 * omega^2 
        + eta2 / (2 * NUCLEON_MASS^2) * gsig^2 * mw^2 * sigma * omega^2
        + etarho / (2 * NUCLEON_MASS) * gsig * mrho^2 * rho^2 + eta1rho / (2 * NUCLEON_MASS^2) * gsig^2 * mrho^2 * sigma * rho^2)]
    push!(val, -gw * (_density(kfp) + _density(kfn)) + mw^2*omega + eta1/NUCLEON_MASS*gsig*mw^2*sigma*omega + eta2/(2*NUCLEON_MASS^2)*gsig^2*mw^2*sigma^2*omega
        + eta2rho/(2*NUCLEON_MASS^2)*gw^2*mrho^2*omega*rho^2 + 1/6*zeta0*gw^2*omega^3)
    push!(val, -1/2 * grho * (_density(kfp) - _density(kfn)) + mrho^2 * rho + etarho / NUCLEON_MASS * gsig * mrho^2 * sigma*rho 
        + eta1rho / (2 * NUCLEON_MASS^2) * gsig^2 * mrho^2 * sigma^2*rho + eta2rho / (2 * NUCLEON_MASS^2) * gw^2 * mrho^2 * omega^2 * rho)
    #apply beta equilibrium
    push!(val, sqrt(kfn^2 + mn_star^2) + gw * omega - 1/2 * grho * rho - mub)
    push!(val, sqrt(kfp^2 + mp_star^2) + gw * omega + 1/2 * grho * rho - mub + mue)
    #apply charge neutrality
    if mue > MUON_MASS
        push!(val, mue^3 + (mue^2 - MUON_MASS^2)^(3/2) - kfp^3)
    else
        push!(val, mue - kfp)
    end
    val
end

function _rmf_energy_density(fields, constants, no_leptons = false) #energy_density in RMF nuclear phase, includes leptons
    msig, mw, mrho, gsig, gw, grho, kappa3, kappa4, zeta0, eta1, eta2, etarho, eta1rho, eta2rho = constants
    sigma, omega, rho, kfp, kfn, mue = fields

    #effective nucleon masses
    mp_star = PROTON_MASS - gsig * sigma
    mn_star = NEUTRON_MASS - gsig * sigma

    #energy density of nucleons and mesons
    nucleon_energy_density = _fermion_energy_density(kfn, mn_star) + _fermion_energy_density(kfp, mp_star)
    meson_energy_density = (1 / 2 * (msig^2 * sigma^2 + mw^2 * omega^2 + mrho^2 * rho^2) 
        + kappa3 / (6 * NUCLEON_MASS) * gsig * msig^2 * sigma^3 + kappa4 / (24 * NUCLEON_MASS^2) * gsig^2 * msig^2 * sigma^4
        + zeta0 / 8 * gw^2 * omega^4 + eta1 / (2 * NUCLEON_MASS) * gsig * mw^2 * sigma * omega^2 + eta2 / (4 * NUCLEON_MASS^2) * gsig^2 * mw^2 * sigma^2 * omega^2
        + etarho / (2 * NUCLEON_MASS) * gsig * mrho^2 * sigma * rho^2 + eta1rho / (4 * NUCLEON_MASS^2) * gsig^2 * mrho^2 * sigma^2 * rho^2 
        + 3 * eta2rho / (4 * NUCLEON_MASS^2) * gw^2 * mrho^2 * omega^2 * rho^2)
    
    #add energy density of leptons
    if no_leptons
        nucleon_energy_density + meson_energy_density
    else 
        nucleon_energy_density + meson_energy_density + _lepton_energy_density(mue)
    end
end
    
function _rmf_pressure(fields, constants, no_leptons = false) #pressure in the RMF nuclear phase, includes leptons
    msig, mw, mrho, gsig, gw, grho, kappa3, kappa4, zeta0, eta1, eta2, etarho, eta1rho, eta2rho = constants
    sigma, omega, rho, kfp, kfn, mue = fields

    #effective nucleon masses
    mp_star = PROTON_MASS - sigma * gsig
    mn_star = NEUTRON_MASS - gsig * sigma

    #pressure from nucleons and mesons
    nucleon_pressure = _fermion_pressure(kfp, mp_star) + _fermion_pressure(kfn, mn_star)
    meson_pressure = (1 / 2 * (-msig^2 * sigma^2 + mw^2 * omega^2 + mrho^2 * rho^2) 
        - kappa3 / (6 * NUCLEON_MASS) * gsig * msig^2 * sigma^3 - kappa4 / (24 * NUCLEON_MASS^2) * gsig^2 * msig^2 * sigma^4
        + zeta0 / 24 * gw^2 * omega^4 + eta1 / (2 * NUCLEON_MASS) * gsig * mw^2 * sigma * omega^2 + eta2 / (4 * NUCLEON_MASS^2) * gsig^2 * mw^2 * sigma^2 * omega^2
        + etarho / (2 * NUCLEON_MASS) * gsig * mrho^2 * sigma * rho^2 + eta1rho / (4 * NUCLEON_MASS^2) * gsig^2 * mrho^2 * sigma^2 * rho^2 
        + eta2rho / (4 * NUCLEON_MASS^2) * gw^2 * mrho^2 * omega^2 * rho^2)

    #add pressure from leptons
    if no_leptons
        nucleon_pressure + meson_pressure
    else
        nucleon_pressure + meson_pressure + _lepton_pressure(mue)
    end
end

struct RMFModel
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

#calculates eos data given mub
function eos_data_rmf(model::RMFModel, mub, nuc_guess = [0.3, 0.2, -0.05, 0.8, 2, 0.8])
    #nuc_guess is in the form [sigma,omega,rho,kfp,kfn,mub,mue]
    #vector of model to pass to RMF functions
    constants = [model.msig, model.mw, model.mrho, model.gsig, model.gw, model.grho, model.kappa3, model.kappa4, model.zeta0,
        model.eta1, model.eta2, model.etarho, model.eta1rho, model.eta2rho]
        #solve equations of motion and store data
    fields_temp = nlsolve(x -> _rmf_eom_neutral(x, mub, constants), nuc_guess)

    kfp = fields_temp.zero[4]
    kfn = fields_temp.zero[5]
    mue = fields_temp.zero[6]
    #calculate eos data from solutions of equations of motion
    energy_density = _rmf_energy_density(fields_temp.zero, constants)
    pressure = _rmf_pressure(fields_temp.zero, constants)
    mp_star = PROTON_MASS - model.gsig * fields_temp.zero[1]
    mn_star = NEUTRON_MASS - model.gsig * fields_temp.zero[1]
    
    [kfp, kfn, mue, energy_density, pressure, mp_star, mn_star]
end

iufsu_star_constants = RMFModel(NUCLEON_MASS * 0.543, NUCLEON_MASS * 0.8331, NUCLEON_MASS * 0.8198, 4 * pi * 0.8379, 
    4 * pi * 1.0666, 4 * pi * 0.9889, 1.1418, 1.0328, 5.3895, 0, 0, 0, 0, 41.3066, "arXiv:1204.2644")
