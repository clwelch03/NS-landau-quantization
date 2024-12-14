const ALPHA_N = 1.13
const BETA_N = 0.68
const NEUTRON_MASS = 939.57
const PROTON_MASS = 938.27

const ERG_TO_MEV = 624151
const RECIP_CM_TO_MEV = 1.973*10^-11
const RECIP_S_TO_MEV = 6.58*10^-22
const G15_TO_FM = 5.00934e-4
const GAUSS_TO_MEV2 = 1.95e-14
const MEV4_TO_CGS = 1 / (ERG_TO_MEV * RECIP_CM_TO_MEV^3 * RECIP_S_TO_MEV)

const G_A = 1.26 
const G_V = 1
const F_2 = 3.7 # duan/qian form factor
const G_F = 1.17e-11 # Fermi coupling constant, MeV^-2
const G_F_BY = 1.44e-49 # Fermi coupling constant for Baiko and Yakovlev, erg cm^3
const CABIBBO_ANGLE = 13.02 * pi/180 # Cabibbo angle, radians
const COS2_CABIBBO_ANGLE = cos(CABIBBO_ANGLE)^2
const ELEM_CHARGE = sqrt(4.0 * pi / 137)
const NEUTRON_AMM = -1.913*ELEM_CHARGE/(2*NEUTRON_MASS)
const PROTON_AMM = 2.793*ELEM_CHARGE/(2*NEUTRON_MASS)
const PROTON_G_MINUS_2 = 3.586
const ELECTRON_G_MINUS_2 = 0.002323
const ELECTRON_MASS = 0.511
const ANALYTICAL_PREFACTOR = (ELEM_CHARGE*G_F^2*COS2_CABIBBO_ANGLE) * 457.0/1290240 * pi
const NUMERICAL_PREFACTOR = (ELEM_CHARGE*G_F^2*COS2_CABIBBO_ANGLE)/(256*pi^5)
const TEMP_STEPS = 15;