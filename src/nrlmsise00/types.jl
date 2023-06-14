# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Types related to the NRLMSISE-00 atmospheric model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export Nrlmsise00Flags, Nrlmsise00Output

"""
    struct Nrlmsise00Flags

Flags to configure NRLMSISE-00.

# Fields

- `F10_Mean::Bool`: F10.7 effect on mean.
- `time_independent::Bool`: Independent of time.
- `sym_annual::Bool`: Symmetrical annual.
- `sym_semiannual::Bool`: Symmetrical semiannual.
- `asym_annual::Bool`: Asymmetrical annual.
- `asyn_semiannual::Bool`: Asymmetrical semiannual.
- `diurnal::Bool`: Diurnal.
- `semidiurnal::Bool`: Semidiurnal.
- `daily_ap::Bool`: Daily AP.
- `all_ut_long_effects::Bool`: All UT/long effects.
- `longitudinal::Bool`: Longitudinal.
- `ut_mixed_ut_long::Bool`: UT and mixed UT/long.
- `mixed_ap_ut_long::Bool`: Mixed AP/UT/long.
- `terdiurnal::Bool`: Terdiurnal.
- `departures_from_eq::Bool`: Departures from diffusive equilibrium.
- `all_tinf_var::Bool`: All TINF variations.
- `all_tlb_var::Bool`: All TLB variations.
- `all_tn1_var::Bool`: All TN1 variations.
- `all_s_var::Bool`: All S variations.
- `all_tn2_var::Bool`: All TN2 variations.
- `all_nlb_var::Bool`: All NLB variations.
- `all_tn3_var::Bool`: All TN3 variations.
- `turbo_scale_height::Bool`: Turbo scale height variations.
"""
Base.@kwdef struct Nrlmsise00Flags
    F10_Mean::Bool            = true
    time_independent::Bool    = true
    sym_annual::Bool          = true
    sym_semiannual::Bool      = true
    asym_annual::Bool         = true
    asym_semiannual::Bool     = true
    diurnal::Bool             = true
    semidiurnal::Bool         = true
    daily_ap::Bool            = true
    all_ut_long_effects::Bool = true
    longitudinal::Bool        = true
    ut_mixed_ut_long::Bool    = true
    mixed_ap_ut_long::Bool    = true
    terdiurnal::Bool          = true
    departures_from_eq::Bool  = true
    all_tinf_var::Bool        = true
    all_tlb_var::Bool         = true
    all_tn1_var::Bool         = true
    all_s_var::Bool           = true
    all_tn2_var::Bool         = true
    all_nlb_var::Bool         = true
    all_tn3_var::Bool         = true
    turbo_scale_height::Bool  = true
end

"""
    mutable struct Nrlmsise00Structure{T<:Number}

Structure with the configuration parameters for NRLMSISE-00 model.
"""
mutable struct Nrlmsise00Structure{T<:Number}
    # Inputs
    # ======================================================================================

    year::Int
    doy::Int
    sec::T
    h::T
    ϕ_gd::T
    λ::T
    lst::T
    F10ₐ::T
    F10::T
    ap::T
    ap_array::Vector{T}
    use_ap_array::Bool
    flags::Nrlmsise00Flags

    # Auxiliary variables to improve code performance
    # ======================================================================================

    r_lat::T
    g_lat::T
    df::T
    dfa::T
    plg::Matrix{T}
    ctloc::T
    stloc::T
    c2tloc::T
    s2tloc::T
    c3tloc::T
    s3tloc::T

    # In the original source code, it has 4 components, but only 1 is used.
    apt::T
    apdf::T
    dm28::T

    # The original code declared all the `meso_*` vectors as global variables.  However,
    # only two values really need to be shared between the functions `gts7` and `gtd7`.
    meso_tn1_5::T
    meso_tgn1_2::T
end

"""
    struct Nrlmsise00Output{T<:Number}

Output structure for NRLMSISE00 model.

# Fields

- `total_density::T`: Total mass density [kg / m³].
- `temperature`: Temperature at the selected altitude [K].
- `exospheric_temperature`: Exospheric temperature [K].
- `N_number_density`: Nitrogen number density [1 / m³].
- `N2_number_density`: N₂ number density [1 / m³].
- `O_number_density`: Oxygen number density [1 / m³].
- `aO_number_density`: Anomalous Oxygen number density [1 / m³].
- `O2_number_density`: O₂ number density [1 / m³].
- `H_number_density`: Hydrogen number density [1 / m³].
- `He_number_density`: Helium number density [1 / m³].
- `Ar_number_density`: Argon number density [1 / m³].

# Remarks

Anomalous oxygen is defined as hot atomic oxygen or ionized oxygen that can become
appreciable at high altitudes (`> 500 km`) for some ranges of inputs, thereby affection drag
on satellites and debris. We group these species under the term **Anomalous Oxygen**, since
their individual variations are not presently separable with the drag data used to define
this model component.
"""
struct Nrlmsise00Output{T<:Number}
    total_density::T
    temperature::T
    exospheric_temperature::T
    N_number_density::T
    N2_number_density::T
    O_number_density::T
    aO_number_density::T
    O2_number_density::T
    H_number_density::T
    He_number_density::T
    Ar_number_density::T
end
