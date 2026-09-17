## Description #############################################################################
#
# Constants for the Jacchia-Bowman 2008 Atmospheric Model.
#
############################################################################################

# Constants to compute the `ΔTc` correction.
const _JB2008_B = (
    -0.457512297e+01,
    -0.512114909e+01,
    -0.693003609e+02,
    0.203716701e+03,
    0.703316291e+03,
    -0.194349234e+04,
    0.110651308e+04,
    -0.174378996e+03,
    0.188594601e+04,
    -0.709371517e+04,
    0.922454523e+04,
    -0.384508073e+04,
    -0.645841789e+01,
    0.409703319e+02,
    -0.482006560e+03,
    0.181870931e+04,
    -0.237389204e+04,
    0.996703815e+03,
    0.361416936e+02,
)

const _JB2008_C = (
    -0.155986211e+02,
    -0.512114909e+01,
    -0.693003609e+02,
    0.203716701e+03,
    0.703316291e+03,
    -0.194349234e+04,
    0.110651308e+04,
    -0.220835117e+03,
    0.143256989e+04,
    -0.318481844e+04,
    0.328981513e+04,
    -0.135332119e+04,
    0.199956489e+02,
    -0.127093998e+02,
    0.212825156e+02,
    -0.275555432e+01,
    0.110234982e+02,
    0.148881951e+03,
    -0.751640284e+03,
    0.637876542e+03,
    0.127093998e+02,
    -0.212825156e+02,
    0.275555432e+01,
)

# F(z) global model values, 1997 - 2006 fit.
const _JB2008_FZM = (0.2689e+00, -0.1176e-01, 0.2782e-01, -0.2782e-01, 0.3470e-03)

# G(t) global model values, 1997 - 2006 fit.
const _JB2008_GTM = (
    -0.3633e+00,
    0.8506e-01,
    0.2401e+00,
    -0.1897e+00,
    -0.2554e+00,
    -0.1790e-01,
    0.5650e-03,
    -0.6407e-03,
    -0.3418e-02,
    -0.1252e-02,
)

# Coefficients for high altitude density correction.
const _JB2008_CHT = (0.22e0, -0.20e-02, 0.115e-02, -0.211e-05)

# Altitude bounds of the model [m].
const _JB2008_H_MIN = 90_000
const _JB2008_H_MAX = 3_000_000

# Physical constants and parameters of the model.
const _JB2008_CONSTANTS = (;
    T₁    = 183.0,      # ............................... Temperature at the lower bound [K]
    z₁    = 90.0,       # ................................. Altitude of the lower bound [km]
    zx    = 125.0,      # ............................ Altitude of the inflection point [km]
    Ra    = 6356.766,   # ......................................... Mean Earth radius [km]
    g₀    = 9.80665,    # .................................. Gravity at Earth surface [m/s²]
    Rstar = 8314.32,    # ............... Universal gas constant (mks) [joules / (K . kmol)]
    ρ₁    = 3.46e-6,    # ........................................ Density at `z₁` [kg / m³]
    A     = 6.02257e26, # ..................... Avogadro's constant (mks) [molecules / kmol]

    # == Assumed Sea-Level Composition =====================================================

    Mb₀  = 28.960,
    q₀N₂ = 0.78110,
    q₀O₂ = 0.20955,
    q₀Ar = 9.3400e-3,
    q₀He = 1.2890e-5,

    # == Molecular Weights of Each Specie [kg / kmol] ======================================

    MN₂ = 28.0134,
    MO₂ = 31.9988,
    MO  = 15.9994,
    MAr = 39.9480,
    MHe = 4.0026,
    MH  = 1.00797,

    # == Thermal Diffusion Coefficient of Each Specie ======================================

    α_N₂ = 0.0,
    α_O₂ = 0.0,
    α_O  = 0.0,
    α_Ar = 0.0,
    α_He = -0.38,
    α_H  = 0.0,

    # == Integration Step Sizes ============================================================
    #
    # `R1` is used between 90 km and 500 km, `R2` in the hydrogen integration from the
    # altitude to 500 km when the altitude is lower than 500 km, and `R3` above 500 km.

    R1 = 0.010,
    R2 = 0.025,
    R3 = 0.075,
)
