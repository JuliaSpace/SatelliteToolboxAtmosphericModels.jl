## Description #############################################################################
#
# The Jacchia-Roberts 1971 Atmospheric Model.
#
## References ##############################################################################
#
# [1] Roberts, C. R (1971). An analytic model for upper atmosphere densities based upon
#     Jacchia's 1970 models.
#
# [2] Jacchia, L. G (1970). New static models of the thermosphere and exosphere with
#     empirical temperature profiles. SAO Special Report #313.
#
# [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
# [4] Long, A. C., Cappellari Jr., J. O., Velez, C. E., Fuchs, A. J (editors) (1989).
#     Goddard Trajectory Determination System (GTDS) Mathematical Theory (Revision 1).
#     FDD/552-89/0001 and CSC/TR-89/6001.
#
# [5] General Mission Analysis Tool (GMAT) R2018a source code.
#     https://sourceforge.net/projects/gmat/
#
############################################################################################

export jr1971

"""
    jr1971(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number) -> JR1971Output
    jr1971(
        instant::DateTime,
        ϕ_gd::Number,
        λ::Number,
        h::Number,
        F10::Number,
        F10ₐ::Number,
        Kp::Number,
    ) -> JR1971Output
    jr1971(jd::Number, ϕ_gd::Number, λ::Number, h::Number) -> JR1971Output
    jr1971(
        jd::Number,
        ϕ_gd::Number,
        λ::Number,
        h::Number,
        F10::Number,
        F10ₐ::Number,
        Kp::Number,
    ) -> JR1971Output

Compute the atmospheric density using the Jacchia-Roberts 1971 model.

If we omit all space indices, the system tries to obtain them automatically for the selected
day `jd` or `instant`. However, the indices must be already initialized using the function
`SpaceIndices.init()`.

The function throws an `ArgumentError` if the altitude `h` is outside the interval
[90, 3000] km.

# Arguments

- `jd::Number`: Julian day to compute the model.
- `instant::DateTime`: Instant to compute the model represented using `DateTime`.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `h::Number`: Altitude [m].
- `F10::Number`: 10.7-cm solar flux adjusted to 1 AU [sfu], as used to fit the Jacchia
    models.
- `F10ₐ::Number`: 10.7-cm averaged solar flux adjusted to 1 AU, 81-day centered on input
    time [sfu].
- `Kp::Number`: Kp geomagnetic index with a delay of 3 hours.

# Returns

- `JR1971Output`: Structure containing the results obtained from the model. Its element
    type is the promotion of the types of the numeric inputs.
"""
function jr1971(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number)
    return jr1971(datetime2julian(instant), ϕ_gd, λ, h)
end

function jr1971(jd::Number, ϕ_gd::Number, λ::Number, h::Number)
    _check_altitude(h, _JR1971_H_MIN, _JR1971_H_MAX)

    # Get the data in the desired Julian Day. The Jacchia models were fitted with the
    # 10.7-cm flux adjusted to 1 AU, so we must not use the observed values here.
    F10  = space_index(Val(:F10adj), jd)
    F10ₐ = _f10_81day_mean(Val(:F10adj), jd)

    # For the Kp, we must obtain the index using a 3-hour delay, considering the Kp constant
    # inside each 3-hour interval provided by the space index vector.
    Kp = _kp_3h(julian2datetime(jd) - Hour(3))

    @debug """
    JR1971 - Fetched Space Indices
      Daily F10.7           : $(F10) sfu
      81-day averaged F10.7 : $(F10ₐ) sfu
      3-hour delayed Kp     : $(Kp)
    """

    return jr1971(jd, ϕ_gd, λ, h, F10, F10ₐ, Kp)
end

function jr1971(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number,
    F10::Number,
    F10ₐ::Number,
    Kp::Number,
)
    return jr1971(datetime2julian(instant), ϕ_gd, λ, h, F10, F10ₐ, Kp)
end

function jr1971(
    jd::JT,
    ϕ_gd::PT,
    λ::LT,
    h::HT,
    F10::FT,
    F10ₐ::FT2,
    Kp::KT,
) where {
    JT <: Number,
    PT <: Number,
    LT <: Number,
    HT <: Number,
    FT <: Number,
    FT2 <: Number,
    KT <: Number,
}
    _check_altitude(h, _JR1971_H_MIN, _JR1971_H_MAX)

    RT = float(promote_type(JT, PT, LT, HT, FT, FT2, KT))

    # == Constants =========================================================================

    Rstar = _JR1971_CONSTANTS.Rstar
    Av    = _JR1971_CONSTANTS.Av
    Ra    = _JR1971_CONSTANTS.Ra
    g₀    = _JR1971_CONSTANTS.g₀
    M₀    = _JR1971_CONSTANTS.M₀
    z₁    = _JR1971_CONSTANTS.z₁
    z₂    = _JR1971_CONSTANTS.z₂
    T₁    = _JR1971_CONSTANTS.T₁
    M₁    = _JR1971_CONSTANTS.M₁
    ρ₁    = _JR1971_CONSTANTS.ρ₁
    Tx_fit = _JR1971_CONSTANTS.Tx_fit
    zx    = _JR1971_CONSTANTS.zx
    Mi    = _JR1971_CONSTANTS.Mi
    αi    = _JR1971_CONSTANTS.αi
    μi    = _JR1971_CONSTANTS.μi
    Aa    = _JR1971_CONSTANTS.Aa
    Ca    = _JR1971_CONSTANTS.Ca
    la    = _JR1971_CONSTANTS.la
    ζ     = _JR1971_CONSTANTS.ζ
    δij   = _JR1971_CONSTANTS.δij

    # == Auxiliary variables ===============================================================

    Ra² = Ra * Ra # ........................................ Mean Earth radius squared [km²]

    # == Preliminaries =====================================================================

    # Convert the altitude from [m] to [km].
    h /= 1000

    # Compute the Sun declination, the Sun right ascension, and the right ascension of the
    # selected location [rad].
    δs, Ωs, Ωp = _sun_geometry(jd, λ)

    # Compute the hour angle at the selected location, which is the angle measured at the XY
    # plane between the right ascension of the selected position and the right ascension of
    # the Sun.
    H = Ωp - Ωs

    ########################################################################################
    #                                      Algorithm                                       #
    ########################################################################################

    # == Exospheric Temperature ============================================================

    # -- Diurnal Variation -----------------------------------------------------------------

    # Eq. 14 [2], Section B.1.1 [3]
    #
    # Nighttime minimum of the global exospheric temperature distribution when the planetary
    # geomagnetic index Kp is zero.

    ΔF10 = F10 - F10ₐ
    Tc   = 379 + 3.24F10ₐ + 1.3ΔF10

    # Eq. 15 [2], Section B.1.1 [3]

    η = abs(ϕ_gd - δs) / 2
    θ = abs(ϕ_gd + δs) / 2

    # Eq. 16 [2], Section B.1.1 [3]

    τ = rem2pi(H + deg2rad(-37 + 6sin(H + deg2rad(43))), RoundNearest)

    # Eq. 17 [2], Section B.1.1 [3]

    C  = cos(η)^2.2
    S  = sin(θ)^2.2
    Tl = Tc * (1 + 0.3 * (S + (C - S) * cos(τ / 2)^3))

    # == Variations with Geomagnetic Activity ==============================================

    # Eq. 18 or Eq. 20 [2], Section B.1.1 [3]

    ΔT∞ = (h < 200) ? 14Kp + 0.02exp(Kp) : 28Kp + 0.03exp(Kp)

    # -- Section B.1.1 [3] -----------------------------------------------------------------
    #
    # Compute the local exospheric temperature with the geomagnetic storm effect.

    T∞ = Tl + ΔT∞

    # == Temperature at the Desired Altitude ===============================================

    # -- Section B.1.1 [3] -----------------------------------------------------------------
    #
    # Compute the temperature at inflection point `zx`.
    #
    # The values at [1, p. 369] are from an old version of Jacchia 1971 model. We will use
    # the new values available at [2].

    Tx = Tx_fit[1] + Tx_fit[2] * T∞ + Tx_fit[3] * exp(Tx_fit[4] * T∞)

    # Compute the temperature at desired point.
    Tz = _jr1971_temperature(h, Tx, T∞)

    # == Corrections to the Density (Eqs. 4-96 to 4-101 [3]) ===============================

    # -- Geomagnetic Effect, Eq. B-7 [3] ---------------------------------------------------

    Δlog₁₀ρ_g = h < 200 ? RT(0.012Kp + 1.2e-5exp(Kp)) : zero(RT)

    # -- Semi-annual Variation, Section B.1.3 [3] ------------------------------------------

    # Number of tropical years since January 1, 1958.
    Φ = (jd - 2436204.5) / 365.2422

    τ_sa = Φ + 0.09544 * ((1 / 2 * (1 + sin(2π * Φ + 6.035)))^(1.65) - 1 / 2)
    f_z  = (5.876e-7h^2.331 + 0.06328) * exp(-0.002868h)
    g_t  = 0.02835 + (0.3817 + 0.17829sin(2π * τ_sa + 4.137)) * sin(4π * τ_sa + 4.259)

    Δlog₁₀ρ_sa = f_z * g_t

    # -- Seasonal Latitudinal Variation, Section B.1.3 [3] ---------------------------------

    sin_ϕ_gd = sin(ϕ_gd)
    abs_sin_ϕ_gd = abs(sin_ϕ_gd)

    Δlog₁₀ρ_lt =
        0.014 *
        (h - 90) *
        exp(-0.0013 * (h - 90)^2) *
        sin(2π * Φ + 1.72) *
        sin_ϕ_gd *
        abs_sin_ϕ_gd

    # -- Total Correction, Eq. B-10 [4] ----------------------------------------------------

    Δlog₁₀ρ_c = Δlog₁₀ρ_g + Δlog₁₀ρ_lt + Δlog₁₀ρ_sa
    Δρ_c = exp10(Δlog₁₀ρ_c)

    # == Density ===========================================================================

    if h <= z₂

        # == Altitudes Between 90 km and 100 km ============================================
        #
        # The closed-form solution of the barometric equation in this region presented in
        # [1] uses partial fractions of the temperature polynomial. Its assembly, as also
        # implemented in [5], leads to densities almost constant between 90 km and 100 km
        # and a discontinuity of a factor of about 6 at 100 km. Hence, we integrate the
        # barometric equation numerically using the Gauss-Legendre quadrature, which is
        # exact to the round-off for the smooth integrand in this short interval.

        # Integral of g M / (R T) between `z₁` and `h` using the 8-point Gauss-Legendre
        # quadrature. Notice that `g` [m / s²], `M` [g / mol], `R` [J / (K . mol)], and the
        # altitude [km] lead to a dimensionless integral.
        Δz = (h - z₁) / 2
        zm = (h + z₁) / 2
        int = zero(RT)

        for i in 1:8
            z   = zm + Δz * _GAUSS_LEGENDRE_8_NODES[i]
            g   = g₀ * Ra² / (Ra + z)^2
            M   = _jr1971_mean_molecular_mass(z)
            T   = _jr1971_temperature(z, Tx, T∞)
            int += _GAUSS_LEGENDRE_8_WEIGHTS[i] * g * M / T
        end

        int *= Δz / Rstar

        # -- Compute the Density, eq. 13 [1] -----------------------------------------------

        Mz = _jr1971_mean_molecular_mass(h)
        ρ  = ρ₁ * Δρ_c * Mz * T₁ / (M₁ * Tz) * exp(-int)

        # The composition is considered constant in this region. The mass density of each
        # species is obtained from the constituent fractions `μi` as in the region between
        # 100 km and 125 km, so that the species mass densities sum to the total density.
        return _jr1971_output(
            RT,
            Tz,
            T∞,
            ρ * μi.N₂ * Mi.N₂ / M₀,
            ρ * μi.O₂ * Mi.O₂ / M₀,
            ρ * μi.O * Mi.O / M₀,
            ρ * μi.Ar * Mi.Ar / M₀,
            ρ * μi.He * Mi.He / M₀,
            ρ * μi.H * Mi.H / M₀,
        )

    elseif h <= zx

        # == Altitudes Between 100 km and 125 km ===========================================

        # First, we need to find the roots of the polynomial:
        #
        #   P(Z) = c₀ + c₁ ⋅ z + c₂ ⋅ z² + c₃ ⋅ z³ + z⁴
        c₀ = (35^4 * Tx / (Tx - T₁) + Ca[1]) / Ca[5]
        c₁ = Ca[2] / Ca[5]
        c₂ = Ca[3] / Ca[5]
        c₃ = Ca[4] / Ca[5]

        r₁, r₂, x, y = _jr1971_roots(c₀, c₁, c₂, c₃)

        # -- k, [1. p. 371] ----------------------------------------------------------------

        k = -g₀ / (Rstar * (Tx - T₁))

        # -- U(ν), V(ν), W(ν), and X — see _jr1971_U, _jr1971_V, _jr1971_W -----------------

        x²_plus_y² = x * x + y * y
        X = -2r₁ * r₂ * Ra * (Ra² + 2x * Ra + x²_plus_y²)

        # First, we need to compute the temperature and density at 100 km.
        T₁₀₀ = _jr1971_temperature(z₂, Tx, T∞)

        # References [3,4] suggest to compute the density using a polynomial fit, so
        # that the computational burden can be reduced:
        #
        ρ₁₀₀ = @evalpoly(T∞, ζ[1], ζ[2], ζ[3], ζ[4], ζ[5], ζ[6], ζ[7]) * M₀

        # Notice that this fit is not exactly consistent with the numerical integration
        # of the barometric equation below 100 km, leading to a small discontinuity (up to
        # about 0.1 %) at 100 km, which is also present in [5].

        # Apply the density correction to the density at 100 km.
        ρ₁₀₀ *= Δρ_c

        # -- Auxiliary Variables, [1, p. 374] ----------------------------------------------

        q₂ = 1 / _jr1971_U(r₁, Ra, x, y, r₁, r₂)
        q₃ = -1 / _jr1971_U(r₂, Ra, x, y, r₁, r₂)
        q₅ = 1 / _jr1971_V(-Ra, x, y, r₁, r₂)
        q₄ =
            (
                1 +
                r₁ * r₂ * (Ra² - x²_plus_y²) * q₅ +
                _jr1971_W(r₁, Ra, x, y, r₁, r₂) * q₂ +
                _jr1971_W(r₂, Ra, x, y, r₁, r₂) * q₃
            ) / X
        q₆ = -q₅ - 2 * (x + Ra) * q₄ - (r₂ + Ra) * q₃ - (r₁ + Ra) * q₂
        q₁ = -2q₄ - q₃ - q₂

        # -- F₃ and F₄, [1, p. 374] --------------------------------------------------------

        log_F₃ =
            q₁ * log((h + Ra) / (z₂ + Ra)) +
            q₂ * log((h - r₁) / (z₂ - r₁)) +
            q₃ * log((h - r₂) / (z₂ - r₂)) +
            q₄ * log((h^2 - 2x * h + x²_plus_y²) / (z₂^2 - 2x * z₂ + x²_plus_y²))

        F₄ =
            q₅ * (h - z₂) / ((h + Ra) * (Ra + z₂)) +
            q₆ / y * atan(y * (h - z₂) / (y^2 + (h - x) * (z₂ - x)))

        # -- Compute the Density of Each Specie [3] ----------------------------------------

        # `f` is defined in [1, p. 371].
        f = 35^4 * Ra² / Ca[5]

        # Notice that we write the power of the temperature ratio as an exponential to
        # evaluate a single transcendental function per species.
        expk = k * f * (log_F₃ + F₄)
        lnT  = log(T₁₀₀ / Tz)
        ρN₂  = ρ₁₀₀ * Mi.N₂ / M₀ * μi.N₂ * exp((1 + αi.N₂) * lnT + Mi.N₂ * expk)
        ρO₂  = ρ₁₀₀ * Mi.O₂ / M₀ * μi.O₂ * exp((1 + αi.O₂) * lnT + Mi.O₂ * expk)
        ρO   = ρ₁₀₀ * Mi.O / M₀ * μi.O * exp((1 + αi.O) * lnT + Mi.O * expk)
        ρAr  = ρ₁₀₀ * Mi.Ar / M₀ * μi.Ar * exp((1 + αi.Ar) * lnT + Mi.Ar * expk)
        ρHe  = ρ₁₀₀ * Mi.He / M₀ * μi.He * exp((1 + αi.He) * lnT + Mi.He * expk)

        return _jr1971_output(RT, Tz, T∞, ρN₂, ρO₂, ρO, ρAr, ρHe, zero(RT))

    else

        # == Altitudes Higher than 125 km ==================================================

        # First, we need to compute the density at 125 km with the corrections.

        # References [3,4] suggest to compute the density using a polynomial fit, so that
        # the computational burden can be reduced:

        ρ₁₂₅_N₂ = Δρ_c * Mi.N₂ * exp10(@evalpoly(T∞, δij.N₂...)) / Av
        ρ₁₂₅_O₂ = Δρ_c * Mi.O₂ * exp10(@evalpoly(T∞, δij.O₂...)) / Av
        ρ₁₂₅_O  = Δρ_c * Mi.O * exp10(@evalpoly(T∞, δij.O...)) / Av
        ρ₁₂₅_Ar = Δρ_c * Mi.Ar * exp10(@evalpoly(T∞, δij.Ar...)) / Av
        ρ₁₂₅_He = Δρ_c * Mi.He * exp10(@evalpoly(T∞, δij.He...)) / Av

        # Notice that this fit leads to a very small discontinuity (about 0.002 %) at
        # 125 km, which is also present in [5].

        # -- Compute `l` According to eq. 4-136 [3] ----------------------------------------

        l = @evalpoly(T∞, la[1], la[2], la[3], la[4], la[5])

        # -- Eq. 25' [1] -------------------------------------------------------------------

        γ   = (g₀ * Ra² / (Rstar * l * T∞) * (T∞ - Tx) / (Tx - T₁) * (zx - z₁) / (Ra + zx))
        γN₂ = γ * Mi.N₂
        γO₂ = γ * Mi.O₂
        γO  = γ * Mi.O
        γAr = γ * Mi.Ar
        γHe = γ * Mi.He

        # -- Eq. 25 [1] --------------------------------------------------------------------

        # Notice that we write the powers of the temperature ratios as exponentials of the
        # logarithms, which are shared by all species, to evaluate a single transcendental
        # function per species.
        lnTx = log(Tx / Tz)
        lnT∞ = log((T∞ - Tz) / (T∞ - Tx))

        ρN₂ = ρ₁₂₅_N₂ * exp((1 + αi.N₂ + γN₂) * lnTx + γN₂ * lnT∞)
        ρO₂ = ρ₁₂₅_O₂ * exp((1 + αi.O₂ + γO₂) * lnTx + γO₂ * lnT∞)
        ρO  = ρ₁₂₅_O * exp((1 + αi.O + γO) * lnTx + γO * lnT∞)
        ρAr = ρ₁₂₅_Ar * exp((1 + αi.Ar + γAr) * lnTx + γAr * lnT∞)
        ρHe = ρ₁₂₅_He * exp((1 + αi.He + γHe) * lnTx + γHe * lnT∞)

        # -- Correction of Seasonal Variations of Helium by Latitude, Eq. 4-101 [3] --------

        Δlog₁₀ρ_He = _jr1971_helium_seasonal_correction(ϕ_gd, δs)

        ρHe *= exp10(Δlog₁₀ρ_He)

        # -- For Altitude Higher than 500 km, We Must Account for H ------------------------

        ρH = zero(ρHe)

        if h > 500
            # Compute the temperature and the H density at 500 km.
            T₅₀₀       = _jr1971_temperature(RT(500), Tx, T∞)
            log₁₀_T₅₀₀ = log10(T₅₀₀)
            ρ₅₀₀_H     = Mi.H / Av * exp10(73.13 - (39.4 - 5.5log₁₀_T₅₀₀) * log₁₀_T₅₀₀)

            # Compute the H density at desired altitude.
            γH = Mi.H * γ
            ρH = Δρ_c * ρ₅₀₀_H * (T₅₀₀ / Tz)^(1 + γH) * ((T∞ - Tz) / (T∞ - T₅₀₀))^γH
        end

        return _jr1971_output(RT, Tz, T∞, ρN₂, ρO₂, ρO, ρAr, ρHe, ρH)
    end
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# The partial-fraction helper functions below are defined at module scope as @inline to
# avoid closures, which incur a runtime allocation from `jl_has_free_typevars` on Julia
# 1.12+.

"""
    _jr1971_U(ν::Number, Ra::Number, x::Number, y::Number, r₁::Number, r₂::Number) -> Number

Compute the partial-fraction helper function `U(ν)` [1, p. 372] given the mean Earth
radius `Ra` [km], the real roots `r₁` and `r₂` [km], and the real and imaginary parts `x`
and `y` [km] of the complex root of the quartic polynomial.
"""
@inline _jr1971_U(ν, Ra, x, y, r₁, r₂) = (ν + Ra)^2 * (ν^2 - 2x * ν + x^2 + y^2) * (r₁ - r₂)

"""
    _jr1971_V(ν::Number, x::Number, y::Number, r₁::Number, r₂::Number) -> Number

Compute the partial-fraction helper function `V(ν)` [1, p. 372] given the real roots `r₁`
and `r₂` [km], and the real and imaginary parts `x` and `y` [km] of the complex root of
the quartic polynomial.
"""
@inline _jr1971_V(ν, x, y, r₁, r₂) = (ν^2 - 2x * ν + x^2 + y^2) * (ν - r₁) * (ν - r₂)

"""
    _jr1971_W(ν::Number, Ra::Number, x::Number, y::Number, r₁::Number, r₂::Number) -> Number

Compute the partial-fraction helper function `W(ν)` [1, p. 372] given the mean Earth
radius `Ra` [km], the real roots `r₁` and `r₂` [km], and the real and imaginary parts `x`
and `y` [km] of the complex root of the quartic polynomial. Notice that the equation for
`W` in [1] is incorrect and the corrected form from [3, 4] is used.
"""
@inline _jr1971_W(ν, Ra, x, y, r₁, r₂) = r₁ * r₂ * Ra * (Ra + ν) * (Ra + (x^2 + y^2) / ν)

############################################################################################

"""
    _jr1971_helium_seasonal_correction(ϕ_gd::Number, δs::Number) -> Number

Compute the base-10 logarithm of the correction of the seasonal variation of the helium
density by latitude (eq. 4-101 [3]) given the geodetic latitude `ϕ_gd` [rad] and the Sun
declination `δs` [rad].
"""
function _jr1971_helium_seasonal_correction(ϕ_gd::Number, δs::Number)
    # Notice that `sign(δs)` is 0 when `δs` is 0 (equinox), which makes the entire
    # expression 0 due to the `abs(δs)` factor. Writing the term as `δs / (2 abs(δs))`, as
    # in [3], would produce a NaN in this case.
    return 0.65 / deg2rad(23.439291) *
           abs(δs) *
           (sin(π / 4 - ϕ_gd * sign(δs) / 2)^3 - 0.35355)
end

"""
    _jr1971_mean_molecular_mass(z::Number) -> Number

Compute the mean molecular mass [g / mol] at the altitude `z` [km] using the empirical
profile in eq. 1 [3, 4], which is valid only between 90 km and 100 km. The caller must
ensure that `z` is inside this range.
"""
function _jr1971_mean_molecular_mass(z::Number)
    Aa = _JR1971_CONSTANTS.Aa
    molecular_mass = @evalpoly(z, Aa[1], Aa[2], Aa[3], Aa[4], Aa[5], Aa[6], Aa[7])

    return molecular_mass
end

"""
    _jr1971_output(
        ::Type{RT},
        Tz::Number,
        T∞::Number,
        ρN₂::Number,
        ρO₂::Number,
        ρO::Number,
        ρAr::Number,
        ρHe::Number,
        ρH::Number
    ) where {RT <: Number} -> JR1971Output{RT}

Build the output structure with element type `RT` given the temperature `Tz` [K] at the
selected altitude, the exospheric temperature `T∞` [K], and the mass densities [g / cm³] of
the species N₂, O₂, O, Ar, He, and H, converting the values to SI units.
"""
function _jr1971_output(
    ::Type{RT},
    Tz::Number,
    T∞::Number,
    ρN₂::Number,
    ρO₂::Number,
    ρO::Number,
    ρAr::Number,
    ρHe::Number,
    ρH::Number,
) where {RT <: Number}
    Av = _JR1971_CONSTANTS.Av
    Mi = _JR1971_CONSTANTS.Mi

    # Total density [g / cm³].
    ρ = ρN₂ + ρO₂ + ρO + ρAr + ρHe + ρH

    # Convert the total density to [kg / m³] and the number densities to [1 / m³].
    return JR1971Output{RT}(
        1000ρ,
        Tz,
        T∞,
        ρN₂ * Av / Mi.N₂ * 1e6,
        ρO₂ * Av / Mi.O₂ * 1e6,
        ρO * Av / Mi.O * 1e6,
        ρAr * Av / Mi.Ar * 1e6,
        ρHe * Av / Mi.He * 1e6,
        ρH * Av / Mi.H * 1e6,
    )
end

"""
    _jr1971_roots(c₀::Number, c₁::Number, c₂::Number, c₃::Number) -> T, T, T, T

Compute the roots of the monic quartic polynomial:

    P(z) = z⁴ + c₃ ⋅ z³ + c₂ ⋅ z² + c₁ ⋅ z + c₀,

which is necessary to compute the density below 125 km. The model theory states that this
polynomial always has two distinct real roots and one complex conjugate pair [1].

The algorithm uses the Ferrari method: the depressed quartic is split into two quadratic
factors whose coefficients are obtained from the largest root of the resolvent cubic,
computed by the Cardano method. The real roots are polished with Newton iterations and the
complex pair is recovered from the Vieta relations, keeping the accuracy close to the
machine precision. Since only closed-form expressions are used, this function does not
allocate and is compatible with automatic differentiation.

# Returns

- `T`: Highest real root `r₁`, where `T` is the promotion of the input types to float.
- `T`: Lowest real root `r₂`.
- `T`: Real part `x` of the complex conjugate pair.
- `T`: Positive imaginary part `y` of the complex conjugate pair.
"""
function _jr1971_roots(c₀::Number, c₁::Number, c₂::Number, c₃::Number)
    c₀, c₁, c₂, c₃ = promote(float(c₀), float(c₁), float(c₂), float(c₃))
    T = typeof(c₀)

    # == Depressed Quartic =================================================================
    #
    # Substituting z = t - c₃ / 4, we obtain the depressed quartic:
    #
    #   t⁴ + p ⋅ t² + q ⋅ t + r = 0 .

    a  = c₃ / 4
    a² = a * a
    p  = c₂ - 6a²
    q  = c₁ - 2c₂ * a + 8a * a²
    r  = c₀ - c₁ * a + c₂ * a² - 3a² * a²

    # == Resolvent Cubic ===================================================================
    #
    # The depressed quartic can be factored as:
    #
    #   (t² + α ⋅ t + β) ⋅ (t² - α ⋅ t + γ) ,
    #
    # in which u = α² is a root of the resolvent cubic:
    #
    #   u³ + 2p ⋅ u² + (p² - 4r) ⋅ u - q² = 0 .
    #
    # Since the cubic is negative at u = 0 and grows unbounded, its largest real root is
    # always non-negative. We compute it using the Cardano method applied to the depressed
    # cubic obtained with u = v - 2p / 3.

    b₂ = 2p
    b₁ = p * p - 4r
    b₀ = -q * q

    P = b₁ - b₂ * b₂ / 3
    Q = (2b₂ * b₂ * b₂ / 9 - b₂ * b₁) / 3 + b₀
    Δ = (Q / 2)^2 + (P / 3)^3

    if Δ >= 0
        # One real root.
        sqrt_Δ = √Δ
        u = cbrt(-Q / 2 + sqrt_Δ) + cbrt(-Q / 2 - sqrt_Δ) - b₂ / 3
    else
        # Three real roots. We take the largest one, which is obtained with k = 0 in the
        # trigonometric solution of the depressed cubic.
        m = √(-P / 3)
        θ = acos(clamp(3Q / (2P * m), -1, 1))
        u = 2m * cos(θ / 3) - b₂ / 3
    end

    u = max(u, zero(T))
    α = √u

    # == Quadratic Factors =================================================================

    if α > √eps(T)
        β = (p + u - q / α) / 2
        γ = (p + u + q / α) / 2

        # One factor contains the two real roots and the other contains the complex
        # conjugate pair. We select them based on the discriminants.
        Δ₁ = u - 4β
        Δ₂ = u - 4γ

        if Δ₁ >= Δ₂
            # Real roots come from t² + α ⋅ t + β = 0.
            sqrt_Δ₁ = √max(Δ₁, zero(T))
            t₊ = (-α + sqrt_Δ₁) / 2
            t₋ = (-α - sqrt_Δ₁) / 2
        else
            # Real roots come from t² - α ⋅ t + γ = 0.
            sqrt_Δ₂ = √max(Δ₂, zero(T))
            t₊ = (α + sqrt_Δ₂) / 2
            t₋ = (α - sqrt_Δ₂) / 2
        end
    else
        # Biquadratic case (q ≈ 0): t² = (-p ± √(p² - 4r)) / 2.
        d  = √max(p * p - 4r, zero(T))
        s₊ = (-p + d) / 2
        t₊ = √max(s₊, zero(T))
        t₋ = -t₊
    end

    # Undo the substitution to obtain the real roots of the original quartic.
    r₁ = t₊ - a
    r₂ = t₋ - a

    # == Newton Polishing ==================================================================

    for _ in 1:2
        f₁  = @evalpoly(r₁, c₀, c₁, c₂, c₃, one(T))
        f₁′ = @evalpoly(r₁, c₁, 2c₂, 3c₃, 4one(T))
        r₁  -= f₁ / f₁′

        f₂  = @evalpoly(r₂, c₀, c₁, c₂, c₃, one(T))
        f₂′ = @evalpoly(r₂, c₁, 2c₂, 3c₃, 4one(T))
        r₂  -= f₂ / f₂′
    end

    if r₁ < r₂
        r₁, r₂ = r₂, r₁
    end

    # == Complex Conjugate Pair ============================================================
    #
    # Using the Vieta relations for the polished real roots, we have:
    #
    #   r₁ + r₂ + 2x = -c₃    and    r₁ ⋅ r₂ ⋅ (x² + y²) = c₀ .

    x = -(c₃ + r₁ + r₂) / 2
    x²_plus_y² = c₀ / (r₁ * r₂)
    y = √(max(x²_plus_y² - x * x, zero(T)))

    return r₁, r₂, x, y
end

"""
    _jr1971_temperature(z::Number, Tx::Number, T∞::Number) -> Number

Compute the temperature [K] at height `z` [km] according to the theory of the model
Jacchia-Roberts 1971 [1, 3, 4] given the temperature `Tx` [K] at the inflection point and
the exospheric temperature `T∞` [K]. The inflection point is considered to be
`z = 125 km`.

The function throws an `ArgumentError` if `z` is lower than 90 km or if `T∞` is negative.
"""
function _jr1971_temperature(z::Number, Tx::Number, T∞::Number)
    T₁ = _JR1971_CONSTANTS.T₁
    z₁ = _JR1971_CONSTANTS.z₁
    zx = _JR1971_CONSTANTS.zx

    # == Check the Parameters ==============================================================

    (z < z₁) && throw(ArgumentError("The altitude must not be lower than $(z₁) km."))
    (T∞ < 0) && throw(ArgumentError("The exospheric temperature must be positive."))

    # == Compute the Temperature at Desired Altitude =======================================

    if z <= zx
        Ca = _JR1971_CONSTANTS.Ca

        aux = @evalpoly(z, Ca[1], Ca[2], Ca[3], Ca[4], Ca[5])
        T   = Tx + (Tx - T₁) / 35^4 * aux
    else
        Ra = _JR1971_CONSTANTS.Ra
        la = _JR1971_CONSTANTS.la

        l = @evalpoly(T∞, la[1], la[2], la[3], la[4], la[5])

        T =
            T∞ -
            (T∞ - Tx) *
            exp(-l * ((Tx - T₁) / (T∞ - Tx)) * ((z - zx) / (zx - z₁)) / (Ra + z))
    end

    return T
end
