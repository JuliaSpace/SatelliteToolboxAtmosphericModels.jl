## Description #############################################################################
#
# Types related to the NRLMSISE-00 atmospheric model.
#
############################################################################################

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
- `asym_semiannual::Bool`: Asymmetrical semiannual.
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
    struct Nrlmsise00Structure{T <: Number, T_AP <: Union{Number, AbstractVector}, T_P <: AbstractMatrix{T}}

Structure with the configuration parameters for NRLMSISE-00 model. `T` is the
floating-number type, `T_AP` is the type of the AP information, which can be a `Number` or
`AbstractVector`, and `T_P` is the type of the matrix with the Legendre associated
functions.

# Fields

- `doy::T`: Day of the year [-].
- `sec::T`: Seconds since the beginning of the day [s].
- `h::T`: Altitude [km].
- `ϕ_gd::T`: Geodetic latitude [°].
- `λ::T`: Longitude [°].
- `lst::T`: Local apparent solar time [h].
- `ap::T_AP`: Magnetic index (daily value or vector with the 3-hour history).
- `flags::Nrlmsise00Flags`: Flags to configure the model.
- `r_lat::T`: Effective Earth radius at the latitude `ϕ_gd` [km].
- `g_lat::T`: Gravity at the latitude `ϕ_gd` [cm / s²].
- `df::T`: Difference between the daily and the 81-day averaged F10.7 flux [sfu].
- `dfa::T`: Difference between the 81-day averaged F10.7 flux and 150 [sfu].
- `plg::T_P`: Unnormalized associated Legendre functions up to degree 7 and order 3, where
    the element `[n + 1, m + 1]` holds the function of degree `n` and order `m` [-].
- `ctloc::T`, `stloc::T`, `c2tloc::T`, `s2tloc::T`, `c3tloc::T`, `s3tloc::T`: Cosine and
    sine of the local solar time and its multiples [-].
- `apt::T`, `apdf::T`: Auxiliary variables of the magnetic activity computed by
    `_globe7` and used by `_glob7s` [-].
- `dm28::T`: N₂ mixed density [1 / cm³] computed by `_gts7` and used by `_gtd7`.
- `meso_tn1_5::T`, `meso_tgn1_2::T`: Temperature [K] and temperature gradient [K / km] at
    the mesopause nodes shared by `_gts7` and `_gtd7`.
"""
struct Nrlmsise00Structure{
    T <: Number,
    T_AP <: Union{Number, AbstractVector},
    T_P <: AbstractMatrix{T},
}
    # == Inputs ============================================================================

    doy::T
    sec::T
    h::T
    ϕ_gd::T
    λ::T
    lst::T
    ap::T_AP
    flags::Nrlmsise00Flags

    # == Auxiliary Variables to Improve Code Performance ===================================

    r_lat::T
    g_lat::T
    df::T
    dfa::T
    plg::T_P
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
    struct Nrlmsise00Output{T <: Number} <: AbstractAtmosphericModelOutput

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
appreciable at high altitudes (`> 500 km`) for some ranges of inputs, thereby affecting drag
on satellites and debris. We group these species under the term **Anomalous Oxygen**, since
their individual variations are not presently separable with the drag data used to define
this model component.
"""
struct Nrlmsise00Output{T <: Number} <: AbstractAtmosphericModelOutput
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

_model_name(::Type{<:Nrlmsise00Output}) = "NRLMSISE-00"
_model_description(::Type{<:Nrlmsise00Output}) = "NRLMSISE-00"

function _show_fields(::Type{<:Nrlmsise00Output})
    return (
        ("Total density", :total_density, _SHOW_FORMAT_DENSITY, "kg / m³"),
        ("Temperature", :temperature, _SHOW_FORMAT_TEMPERATURE, "K"),
        ("Exospheric Temp.", :exospheric_temperature, _SHOW_FORMAT_TEMPERATURE, "K"),
        ("N  number density", :N_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("N₂ number density", :N2_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("O  number density", :O_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("Anomalous O num. den.", :aO_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("O₂ number density", :O2_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("Ar number density", :Ar_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("He number density", :He_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("H  number density", :H_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
    )
end
