# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   The Jacchia-Roberts 1971 Atmospheric Model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Roberts, C. R (1971). An analytic model for upper atmosphere densities based upon
#       Jacchia's 1970 models.
#
#   [2] Jacchia, L. G (1970). New static models of the thermosphere and exosphere with
#       empirical temperature profiles. SAO Special Report #313.
#
#   [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [4] Long, A. C., Cappellari Jr., J. O., Velez, C. E., Fuchs, A. J (editors) (1989).
#       Goddard Trajectory Determination System (GTDS) Mathematical Theory (Revision 1).
#       FDD/552-89/0001 and CSC/TR-89/6001.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export jr1971

"""
    jr1971(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number[, F10::Number, F10ₐ::Number, Kp::Number]) -> JR1971Output{Float64}
    jr1971(jd::Number, ϕ_gd::Number, λ::Number, h::Number[, F10::Number, F10ₐ::Number, Kp::Number]) -> JR1971Output{Float64}

Compute the atmospheric density using the Jacchia-Roberts 1971 model.

If we omit all space indices, the system tries to obtain them automatically for the selected
day `jd` or `instant`. However, the indices must be already initialized using the function
`SpaceIndices.init()`.

# Arguments

- `jd::Number`: Julian day to compute the model.
- `instant::DateTime`: Instant to compute the model represent using `DateTime`.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `h::Number`: Altitude [m].
- `F10::Number`: 10.7-cm solar flux [sfu].
- `F10ₐ::Number`: 10.7-cm averaged solar flux, 81-day centered on input time [sfu].
- `Kp::Number`: Kp geomagnetic index with a delay of 3 hours.

# Returns

- `JR1971Output{Float64}`: Structure containing the results obtained from the model.
"""
function jr1971(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number)
    return jr1971(datetime2julian(instant), ϕ_gd, λ, h)
end

function jr1971(jd::Number, ϕ_gd::Number, λ::Number, h::Number)
    # Get the data in the desired Julian Day.
    F10  = space_index(Val(:F10obs), jd)
    F10ₐ = sum((space_index.(Val(:F10obs), k) for k in (jd - 40):(jd + 40))) / 81

    # For the Kp, we must obtain the index using a 3-hour delay. Thus, we need to obtain the
    # Kp vector first, containing the Kp values for every 3 hours.
    instant = julian2datetime(jd)
    instant_3h_delay = instant - Hour(3)
    Kp_vect = space_index(Val(:Kp), instant_3h_delay)

    # Now, we need to obtain the value for the required instant. In this case, we will
    # consider the Kp constant inside the 3h-interval provided by the space index vector.

    # Get the number of seconds elapsed since the beginning of the day.
    day = Date(instant) |> DateTime
    Δt = Dates.value(instant - day) / 1000

    # Get the index and the Kp value.
    id = clamp(div(Δt, 10_800) + 1, 1, 8)
    Kp = Kp_vect[id]

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
    Kp::Number
)
    return jr1971(datetime2julian(instant), ϕ_gd, λ, h, F10, F10ₐ, Kp)
end

function jr1971(
    jd::Number,
    ϕ_gd::Number,
    λ::Number,
    h::Number,
    F10::Number,
    F10ₐ::Number,
    Kp::Number
)
    # Constants
    # ======================================================================================

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
    zx    = _JR1971_CONSTANTS.zx
    Mi    = _JR1971_CONSTANTS.Mi
    αi    = _JR1971_CONSTANTS.αi
    μi    = _JR1971_CONSTANTS.μi
    Aa    = _JR1971_CONSTANTS.Aa
    Ca    = _JR1971_CONSTANTS.Ca
    la    = _JR1971_CONSTANTS.la
    α     = _JR1971_CONSTANTS.α
    β     = _JR1971_CONSTANTS.β
    ζ     = _JR1971_CONSTANTS.ζ
    δij   = _JR1971_CONSTANTS.δij

    # Auxiliary variables
    # ======================================================================================

    Ra² = Ra * Ra # ........................................ Mean Earth radius squared [km²]

    # Preliminaries
    # ======================================================================================

    # Convert the altitude from [m] to [km].
    h /= 1000

    # Compute the Sun position represented in the inertial reference frame (MOD).
    s_i = sun_position_mod(jd)

    # Compute the Sun declination [rad].
    δs = atan(s_i[3], √(s_i[1] * s_i[1] + s_i[2] * s_i[2]))

    # Compute the Sun right ascension [rad].
    Ωs = atan(s_i[2], s_i[1])

    # Compute the right ascension of the selected location w.r.t. the inertial reference
    # frame.
    Ωp = λ + jd_to_gmst(jd)

    # Compute the hour angle at the selected location, which is the angle measured at the XY
    # plane between the right ascension of the selected position and the right ascension of
    # the Sun.
    H = Ωp - Ωs

    ########################################################################################
    #                                      Algorithm
    ########################################################################################

    #                                Exospheric Temperature
    # ======================================================================================

    # Diurnal variation
    # ======================================================================================

    # Eq. 14 [2], Section B.1.1 [3]
    # --------------------------------------------------------------------------------------
    #
    # Nighttime minimum of the global exospheric temperature distribution when the planetary
    # geomagnetic index Kp is zero.

    ΔF10 = F10 - F10ₐ
    Tc   = 379 + 3.24F10 + 1.3ΔF10

    # Eq. 15 [2], Section B.1.1 [3]
    # --------------------------------------------------------------------------------------

    η = abs(ϕ_gd - δs) / 2
    θ = abs(ϕ_gd + δs) / 2

    # Eq. 16 [2], Section B.1.1 [3]
    # --------------------------------------------------------------------------------------

    τ = H + deg2rad(-37 + 6sin(H + deg2rad(43)))

    # Eq. 17 [2], Section B.1.1 [3]
    # --------------------------------------------------------------------------------------

    C  = cos(η)^2.2
    S  = sin(θ)^2.2
    Tl = Tc * (1 + 0.3 * (S + (C - S) * cos(τ / 2)^3))

    # Variations with geomagnetic activity
    # ======================================================================================

    # Eq. 18 or Eq. 20 [2], Section B.1.1 [3]
    # --------------------------------------------------------------------------------------

    ΔT∞ = (h < 200) ? 14Kp + 0.02exp(Kp) : 28Kp + 0.03exp(Kp)

    # Section B.1.1 [3]
    # --------------------------------------------------------------------------------------
    #
    # Compute the local exospheric temperature with the geomagnetic storm effect.

    T∞ = Tl + ΔT∞

    #                         Temperature at the Desired Altitude
    # ======================================================================================

    # Section B.1.1 [3]
    # --------------------------------------------------------------------------------------
    #
    # Compute the temperature at inflection point `zx`.
    #
    # The values at [1, p. 369] are from an old version of Jacchia 1971 model. We will use
    # the new values available at [2].

    a  =  371.6678
    b  =    0.0518806
    c  = -294.3505
    d  =   -0.00216222
    Tx = a + b * T∞ + c * exp(d * T∞)

    # Compute the temperature at desired point.
    Tz = _jr1971_temperature(h, Tx, T∞)

    #                 Corrections to the Density (Eqs. 4-96 to 4-101 [3])
    # ======================================================================================

    # Geomagnetic effect, Eq. B-7 [3]
    # --------------------------------------------------------------------------------------

    Δlog₁₀ρ_g = h < 200 ? 0.012Kp + 1.2e-5exp(Kp) : 0.0

    # Semi-annual variation, Section B.1.3 [3]
    # --------------------------------------------------------------------------------------

    # Number of days since January 1, 1958.
    Φ = (jd - 2436204.5) / 365.2422

    τ_sa = Φ + 0.09544 * ((1 / 2 * (1 + sin(2π * Φ + 6.035)))^(1.65) - 1 / 2)
    f_z  = (5.876e-7h^2.331 + 0.06328) * exp(-0.002868h)
    g_t  = 0.02835 + (0.3817 + 0.17829sin(2π * τ_sa + 4.137)) * sin(4π * τ_sa + 4.259)

    Δlog₁₀ρ_sa = f_z * g_t

    # Seasonal latitudinal variation, Section B.1.3 [3]
    # --------------------------------------------------------------------------------------

    sin_ϕ_gd = sin(ϕ_gd)
    abs_sin_ϕ_gd = abs(sin_ϕ_gd)

    Δlog₁₀ρ_lt = 0.014 * (h - 90) * exp(-0.0013 * (h - 90)^2) * sin(2π * Φ + 1.72) * sin_ϕ_gd * abs_sin_ϕ_gd

    # Total correction, Eq. B-10 [4]
    # --------------------------------------------------------------------------------------

    Δlog₁₀ρ_c = Δlog₁₀ρ_g + Δlog₁₀ρ_lt + Δlog₁₀ρ_sa
    Δρ_c = 10^(Δlog₁₀ρ_c)

    #                                       Density
    # ======================================================================================

    if h == z₁
        # Compute the total density.
        ρ = ρ₁ * Δρ_c

        # Convert to SI and return.
        return JR1971Output(
            1000ρ,
            Tz,
            T∞,
            (ρ * μi.N₂) * Av / Mi.N₂ * 1e6,
            (ρ * μi.O₂) * Av / Mi.O₂ * 1e6,
            (ρ * μi.O)  * Av / Mi.O  * 1e6,
            (ρ * μi.Ar) * Av / Mi.Ar * 1e6,
            (ρ * μi.He) * Av / Mi.He * 1e6,
            (ρ * μi.H)  * Av / Mi.H  * 1e6
        )

    elseif z₁ < h <= zx

        # First, we need to find the roots of the polynomial:
        #
        #   P(Z) = c₀ + c₁ ⋅ z + c₂ ⋅ z² + c₃ ⋅ z³ + c₄ ⋅ z⁴

        c₀ = (35^4 * Tx / (Tx - T₁) + Ca[1]) / Ca[5]
        c₁ = Ca[2] / Ca[5]
        c₂ = Ca[3] / Ca[5]
        c₃ = Ca[4] / Ca[5]
        c₄ = Ca[5] / Ca[5]
        r₁, r₂, x, y = _jr1971_roots([c₀, c₁, c₂, c₃, c₄])

        # f and k, [1. p. 371]
        # ----------------------------------------------------------------------------------

        f = 35^4 * Ra² / Ca[5]
        k = -g₀ / (Rstar * (Tx - T₁))

        # U(ν), V(ν), and W(ν) functions, [1, p. 372]
        # ----------------------------------------------------------------------------------

        U(ν) = (ν + Ra)^2 * (ν^2 - 2x * ν + x^2 + y^2) * (r₁ - r₂)
        V(ν) = (ν^2 - 2x * ν + x^2 + y^2) * (ν - r₁) * (ν - r₂)

        # The equation for `W(ν)` in [1] was:
        #
        #   W(ν) = r₁ * r₂ * Ra² * (ν + Ra) + (x^2 + y^2) * Ra * (Ra * ν + r₁ * r₂)
        #
        # However, [3,4] mention that this is not correct, and must be replaced by:

        W(ν) = r₁ * r₂ * Ra * (Ra + ν) * (Ra + (x^2 + y^2) / ν)

        # X, [1, p. 372]
        # ----------------------------------------------------------------------------------

        X = -2r₁ * r₂ * Ra * (Ra² + 2x * Ra + x^2 + y^2)

        # The original paper [1] divided into two sections: 90 km to 105 km, and 105 km to
        # 125 km. However, [3] divided into 90 to 100 km, and 100 km to 125 km.

        if h <= z₂

            # Altitudes between 90 km and 100 km
            # ==============================================================================

            # S(z) polynomial, [1, p. 371]
            # ------------------------------------------------------------------------------

            B₀ = α[1] + β[1] * Tx / (Tx - T₁)
            B₁ = α[2] + β[2] * Tx / (Tx - T₁)
            B₂ = α[3] + β[3] * Tx / (Tx - T₁)
            B₃ = α[4] + β[4] * Tx / (Tx - T₁)
            B₄ = α[5] + β[5] * Tx / (Tx - T₁)
            B₅ = α[6] + β[6] * Tx / (Tx - T₁)

            S(z) = @evalpoly(z, B₀, B₁, B₂, B₃, B₄, B₅)

            # Auxiliary variables, [1. p. 372]
            # ------------------------------------------------------------------------------

            p₂ =  S(r₁) / U(r₁)
            p₃ = -S(r₂) / U(r₂)
            p₅ = S(-Ra) / V(-Ra)

            # There is a typo in the fourth term in [1] that was corrected in [3].

            p₄ = (
                B₀ - r₁ * r₂ * Ra² * (B₄ + (2x + r₁ + r₂ - Ra) * B₅) -
                r₁ * r₂ * Ra * (x^2 + y^2) * B₅ + r₁ * r₂ * (Ra² - (x^2 + y^2)) * p₅ +
                W(r₁) * p₂ + W(r₂) * p₃
            ) / X

            p₆ = B₄ + (2x + r₁ + r₂ - Ra ) * B₅ - p₅ - 2(x + Ra) * p₄ - (r₂ + Ra) * p₃ - (r₁ + Ra) * p₂
            p₁ = B₅ - 2p₄ - p₃ - p₂

            # F₁ and F₂, [1, p. 372-373]
            # ------------------------------------------------------------------------------

            log_F₁ = p₁ * log((h + Ra) / (z₁ + Ra)) +
                     p₂ * log((h - r₁) / (z₁ - r₁)) +
                     p₃ * log((h - r₂) / (z₁ - r₂)) +
                     p₄ * log((h^2 - 2x * h  + x^2 + y^2 ) / (z₁^2 - 2x * z₁ + x^2 + y^2))

            # This equation in [4] is wrong, since `f` is multiplying `A₆`. We will use the
            # one in [3].
            F₂ = (h - z₁) * (Aa[7] + p₅ / ((h + Ra) * (z₁ + Ra))) +
                 p₆ / y * atan(y * (h - z₁) / (y^2 + (h - x) * (z₁ - x)))

            # Compute the density, eq. 13 [1]
            # ------------------------------------------------------------------------------

            Mz = _jr1971_mean_molecular_mass(h)
            ρ  = ρ₁ * Δρ_c * Mz * T₁ / (M₁ * Tz) * exp(k * (log_F₁ + F₂))

            # Convert to SI and return.
            return JR1971Output(
                1000ρ,
                Tz,
                T∞,
                (ρ * μi.N₂) * Av / Mi.N₂ * 1e6,
                (ρ * μi.O₂) * Av / Mi.O₂ * 1e6,
                (ρ * μi.O)  * Av / Mi.O  * 1e6,
                (ρ * μi.Ar) * Av / Mi.Ar * 1e6,
                (ρ * μi.He) * Av / Mi.He * 1e6,
                (ρ * μi.H)  * Av / Mi.H  * 1e6,
            )
        else

            # Altitudes between 100 km and 125 km
            # ==============================================================================

            # First, we need to compute the temperature and density at 100 km.
            T₁₀₀ = _jr1971_temperature(z₂, Tx, T∞)

            # References [3,4] suggest to compute the density using a polynomial fit, so
            # that the computational burden can be reduced:
            #
            ρ₁₀₀ = @evalpoly(T∞, ζ[1], ζ[2], ζ[3], ζ[4], ζ[5], ζ[6], ζ[7]) * M₀

            # However, it turns out that this approach leads to discontinuity at 100 km.
            # This was also seen by GMAT [5].
            #
            # TODO: Check if we need to fix this.

            # Apply the density correction to the density at 100 km.
            ρ₁₀₀ *= Δρ_c

            # Auxiliary variables, [1, p. 374]
            # ------------------------------------------------------------------------------

            q₂ =  1 / U(r₁)
            q₃ = -1 / U(r₂)
            q₅ =  1 / V(-Ra)
            q₄ = (1 + r₁ * r₂ * (Ra² - (x^2 + y^2)) * q₅ + W(r₁) * q₂ + W(r₂) * q₃) / X
            q₆ = -q₅ - 2 * (x + Ra) * q₄ - (r₂ + Ra) * q₃ - (r₁ + Ra) * q₂
            q₁ = -2q₄ - q₃ - q₂

            # F₃ and F₄, [1, p. 374]
            # ------------------------------------------------------------------------------

            log_F₃ = q₁ * log((h + Ra) / (z₂ + Ra)) +
                     q₂ * log((h - r₁) / (z₂ - r₁)) +
                     q₃ * log((h - r₂) / (z₂ - r₂)) +
                     q₄ * log((h^2 - 2x*h + x^2 + y^2) / (z₂^2 - 2x*z₂ + x^2 + y^2))

            F₄ = q₅ * (h - z₂) / ((h + Ra) * (Ra + z₂)) +
                 q₆ / y * atan(y * (h - z₂) / (y^2 + (h - x) * (z₂ - x)))

            # Compute the density of each specie [3]
            # ------------------------------------------------------------------------------

            expk = k * f * (log_F₃ + F₄)
            ρN₂  = ρ₁₀₀ * Mi[1] / M₀ * μi[1] * (T₁₀₀ / Tz)^(1 + αi.N₂) * exp(Mi[1] * expk)
            ρO₂  = ρ₁₀₀ * Mi[2] / M₀ * μi[2] * (T₁₀₀ / Tz)^(1 + αi.O₂) * exp(Mi[2] * expk)
            ρO   = ρ₁₀₀ * Mi[3] / M₀ * μi[3] * (T₁₀₀ / Tz)^(1 + αi.O ) * exp(Mi[3] * expk)
            ρAr  = ρ₁₀₀ * Mi[4] / M₀ * μi[4] * (T₁₀₀ / Tz)^(1 + αi.Ar) * exp(Mi[4] * expk)
            ρHe  = ρ₁₀₀ * Mi[5] / M₀ * μi[5] * (T₁₀₀ / Tz)^(1 + αi.He) * exp(Mi[5] * expk)

            # Compute the total density.
            ρ = ρN₂ + ρO₂ + ρO + ρAr + ρHe

            # Convert to SI and return.
            return JR1971Output(
                1000ρ,
                Tz,
                T∞,
                ρN₂ * Av / Mi.N₂ * 1e6,
                ρO₂ * Av / Mi.O₂ * 1e6,
                ρO  * Av / Mi.O  * 1e6,
                ρAr * Av / Mi.Ar * 1e6,
                ρHe * Av / Mi.He * 1e6,
                0.0,
            )
        end

    else

        # Altitudes higher than 125 km
        # ==================================================================================

        # First, we need to compute the density at 125 km with the corrections.

        # References [3,4] suggest to compute the density using a polynomial fit, so that
        # the computational burden can be reduced:

        ρ₁₂₅_N₂  = Δρ_c * Mi[1] * 10^(@evalpoly(T∞, δij.N₂...)) / Av
        ρ₁₂₅_O₂  = Δρ_c * Mi[2] * 10^(@evalpoly(T∞, δij.O₂...)) / Av
        ρ₁₂₅_O   = Δρ_c * Mi[3] * 10^(@evalpoly(T∞, δij.O...))  / Av
        ρ₁₂₅_Ar  = Δρ_c * Mi[4] * 10^(@evalpoly(T∞, δij.Ar...)) / Av
        ρ₁₂₅_He  = Δρ_c * Mi[5] * 10^(@evalpoly(T∞, δij.He...)) / Av

        # However, it turns out that this approach leads to discontinuity at 100 km. This
        # was also seen by GMAT [5].
        #
        # TODO: Check if we need to fix this.

        # Compute `l` according to eq. 4-136 [3]
        # ----------------------------------------------------------------------

        l = @evalpoly(T∞, la[1], la[2], la[3], la[4], la[5])

        # Eq. 25' [1]
        # ----------------------------------------------------------------------

        γ   = (g₀ * Ra² / (Rstar * l * T∞) * (T∞ - Tx) / (Tx - T₁) * (zx - z₁) / (Ra + zx))
        γN₂ = γ * Mi[1]
        γO₂ = γ * Mi[2]
        γO  = γ * Mi[3]
        γAr = γ * Mi[4]
        γHe = γ * Mi[5]

        # Eq. 25 [1]
        # ----------------------------------------------------------------------------------

        ρN₂  = ρ₁₂₅_N₂ * (Tx / Tz)^(1 + αi.N₂ + γN₂) * ((T∞ - Tz) / (T∞ - Tx)) ^ γN₂
        ρO₂  = ρ₁₂₅_O₂ * (Tx / Tz)^(1 + αi.O₂ + γO₂) * ((T∞ - Tz) / (T∞ - Tx)) ^ γO₂
        ρO   = ρ₁₂₅_O  * (Tx / Tz)^(1 + αi.O  + γO ) * ((T∞ - Tz) / (T∞ - Tx)) ^ γO
        ρAr  = ρ₁₂₅_Ar * (Tx / Tz)^(1 + αi.Ar + γAr) * ((T∞ - Tz) / (T∞ - Tx)) ^ γAr
        ρHe  = ρ₁₂₅_He * (Tx / Tz)^(1 + αi.He + γHe) * ((T∞ - Tz) / (T∞ - Tx)) ^ γHe

        # Correction of seasonal variations of helium by latitude, Eq. 4-101 [3]
        # ----------------------------------------------------------------------------------

        Δlog₁₀ρ_He = 0.65 / deg2rad(23.439291) * abs(δs) * (
            sin(π / 4 - ϕ_gd * δs / (2abs(δs)))^3 - 0.35355
        )

        ρHe *= 10^(Δlog₁₀ρ_He)

        # For altitude higher than 500 km, we must account for H
        # ----------------------------------------------------------------------------------

        ρH = 0.0

        if h > 500
            # Compute the temperature and the H density at 500 km.
            T₅₀₀       = _jr1971_temperature(500.0, Tx, T∞)
            log₁₀_T₅₀₀ = log10(T₅₀₀)
            ρ₅₀₀_H     = Mi.H / Av * 10^(73.13 - (39.4 - 5.5log₁₀_T₅₀₀) * log₁₀_T₅₀₀)

            # Compute the H density at desired altitude.
            γH = Mi.H * γ
            ρH = Δρ_c * ρ₅₀₀_H * (T₅₀₀ / Tz)^(1 + γH) * ((T∞ - Tz) / (T∞ - T₅₀₀))^γH
        end

        # Apply the corrections.
        ρ = ρN₂ + ρO₂ + ρO + ρAr + ρHe + ρH

        # Convert to SI and return.
        return JR1971Output(
            1000ρ,
            Tz,
            T∞,
            ρN₂ * Av / Mi.N₂ * 1e6,
            ρO₂ * Av / Mi.O₂ * 1e6,
            ρO  * Av / Mi.O  * 1e6,
            ρAr * Av / Mi.Ar * 1e6,
            ρHe * Av / Mi.He * 1e6,
            ρH  * Av / Mi.H  * 1e6,
        )
    end
end

############################################################################################
#                                    Private Functions
############################################################################################

#   _jr1971_mean_molecular_mass(z::Number) -> Float64
#
# Compute the mean molecular mass at altitude `z` [km] using the empirical profile in eq. 1
# **[3, 4]**.
function _jr1971_mean_molecular_mass(z::Number)
    !(90 <= z <= 100) &&
        @warn("The empirical model for the mean molecular mass is valid only for 90 <= z <= 100 km.")

    Aa = _JR1971_CONSTANTS.Aa
    molecular_mass = @evalpoly(z, Aa[1], Aa[2], Aa[3], Aa[4], Aa[5], Aa[6], Aa[7])

    return molecular_mass
end

#     _jr1971_roots(p::AbstractVector) where T<:Number -> NTuple{4, Float64}
#
# Compute the roots of the polynomial `p` necessary to compute the density below 125 km. It
# returns `r₁`, `r₂`, `x`, and `y`.
function _jr1971_roots(p::AbstractVector)
    # Compute the roots with a first guess.
    r = roots(p, _JR1971_ROOT_GUESS; polish = true)

    # We expect two real roots and two complex roots. Here, we will perform the following
    # processing:
    #
    #   r₁ -> Highest real root.
    #   r₂ -> Lowest real root.
    #   x  -> Real part of the complex root.
    #   y  -> Positive imaginary part of the complex root.

    r₁ = maximum(v -> abs(imag(v)) < 1e-10 ? real(v) : -Inf, r)
    r₂ = minimum(v -> abs(imag(v)) < 1e-10 ? real(v) : +Inf, r)
    c  = findfirst(v -> imag(v) >= 1e-10, r)
    x  = real(r[c])
    y  = abs(imag(r[c]))

    return r₁, r₂, x, y
end

#   _jr1971_temperature(z::Number, Tx::Number, T∞::Number) -> Float64
#
# Compute the temperature [K] at height `z` [km] according to the theory of the model
# Jacchia-Roberts 1971 [1, 3, 4] given the temperature `Tx` [K] at the inflection point and
# the exospheric temperature `T∞` [K].
#
# The inflection point is considered to by `z = 125 km`.
function _jr1971_temperature(z::Number, Tx::Number, T∞::Number)
    T₁ = _JR1971_CONSTANTS.T₁
    z₁ = _JR1971_CONSTANTS.z₁
    zx = _JR1971_CONSTANTS.zx

    # Check the parameters
    # ======================================================================================

    (z  < z₁) && throw(ArgumentError("The altitude must not be lower than $(z₁) km."))
    (T∞ < 0 ) && throw(ArgumentError("The exospheric temperature must be positive."))

    # Compute the temperature at desire altitude
    # ======================================================================================

    if z <= zx
        Ca = _JR1971_CONSTANTS.Ca

        aux = @evalpoly(z, Ca[1], Ca[2], Ca[3], Ca[4], Ca[5])
        T   = Tx + (Tx - T₁) / 35^4 * aux
    else
        Ra = _JR1971_CONSTANTS.Ra
        la = _JR1971_CONSTANTS.la

        l = @evalpoly(T∞, la[1], la[2], la[3], la[4], la[5])

        T = T∞ - (T∞ - Tx) * exp(
            -l * ((Tx - T₁) / (T∞ - Tx)) * ((z - zx) / (zx - z₁)) / (Ra + z)
        )
    end

    return T
end
