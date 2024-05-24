## Description #############################################################################
#
# The Jacchia-Bowman 2008 Atmospheric Model, a product of the Space Environment
# Technologies.
#
# The source code is available at the following git:
#
#     http://sol.spacenvironment.net/jb2008/code.html
#
# For more information about the model, see:
#
#     http://sol.spacenvironment.net/~JB2008/
#
## References ##############################################################################
#
# [1] Bowman, B. R., Tobiska, W. K., Marcos, F. A., Huang, C. Y., Lin, C. S., Burke, W. J
#     (2008). A new empirical thermospheric density model JB2008 using new solar and
#     geomagnetic indices. AIAA/AAS Astrodynamics Specialist Conference, Honolulu, Hawaii.
#
# [2] Bowman, B. R., Tobiska, W. K., Marcos, F. A., Valladares, C (2007). The JB2006
#     empirical thermospheric density model. Journal of Atmospheric and Solar-Terrestrial
#     Physics, v. 70, p. 774-793.
#
# [3] Jacchia, L. G (1970). New static models of the thermosphere and exosphere with
#     empirical temperature profiles. SAO Special Report #313.
#
############################################################################################

export jb2008

"""
    jb2008(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number[, F10::Number, F10ₐ::Number, S10::Number, S10ₐ::Number, M10::Number, M10ₐ::Number, Y10::Number, Y10ₐ::Number, DstΔTc::Number]) -> JB2008Output{Float64}
    jb2008(jd::Number, ϕ_gd::Number, λ::Number, h::Number[, F10::Number, F10ₐ::Number, S10::Number, S10ₐ::Number, M10::Number, M10ₐ::Number, Y10::Number, Y10ₐ::Number, DstΔTc::Number]) -> JB2008Output{Float64}

Compute the atmospheric density using the Jacchia-Bowman 2008 (JB2008) model.

This model is a product of the **Space Environment Technologies**, please, refer to the
following website for more information:

http://sol.spacenvironment.net/JB2008/

If we omit all space indices, the system tries to obtain them automatically for the selected
day `jd` or `instant`. However, the indices must be already initialized using the function
`SpaceIndices.init()`.

# Arguments

- `jd::Number`: Julian day to compute the model.
- `instant::DateTime`: Instant to compute the model represent using `DateTime`.
- `ϕ_gd`: Geodetic latitude [rad].
- `λ`: Longitude [rad].
- `h`: Altitude [m].
- `F10`: 10.7-cm solar flux [sfu] obtained 1 day before `jd`.
- `F10ₐ`: 10.7-cm averaged solar flux using a 81-day window centered on input time obtained
    1 day before `jd`.
- `S10`: EUV index (26-34 nm) scaled to F10.7 obtained 1 day before `jd`.
- `S10ₐ`: EUV 81-day averaged centered index obtained 1 day before `jd`.
- `M10`: MG2 index scaled to F10.7 obtained 2 days before `jd`.
- `M10ₐ`: MG2 81-day averaged centered index obtained 2 day before `jd`.
- `Y10`: Solar X-ray & Ly-α index scaled to F10.7 obtained 5 days before `jd`.
- `Y10ₐ`: Solar X-ray & Ly-α 81-day averaged centered index obtained 5 days before `jd`.
- `DstΔTc`: Temperature variation related to the Dst.

# Returns

- `JB2008Output{Float64}`: Structure containing the results obtained from the model.
"""
function jb2008(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number)
    return jb2008(datetime2julian(instant), ϕ_gd, λ, h)
end

function jb2008(jd::Number, ϕ_gd::Number, λ::Number, h::Number)
    # Get the data in the desired Julian Day considering the tabular time of the model.
    F10    = space_index(Val(:F10obs), jd - 1)
    F10ₐ   = sum((space_index.(Val(:F10obs), k) for k in (jd - 1 - 40):(jd - 1 + 40))) / 81
    S10    = space_index(Val(:S10), jd - 1)
    S10ₐ   = space_index(Val(:S81a), jd - 1)
    M10    = space_index(Val(:M10), jd - 2)
    M10ₐ   = space_index(Val(:M81a), jd - 2)
    Y10    = space_index(Val(:Y10), jd - 5)
    Y10ₐ   = space_index(Val(:Y81a), jd - 5)
    DstΔTc = space_index(Val(:DTC), jd)

    @debug """
    JB2008 - Fetched Space Indices
      Daily F10.7           : $(F10) sfu
      81-day averaged F10.7 : $(F10ₐ) sfu
      Daily S10             : $(S10)
      81-day averaged S10   : $(S10ₐ)
      Daily M10             : $(M10)
      81-day averaged M10   : $(M10ₐ)
      Daily Y10             : $(Y10)
      81-day averaged Y10   : $(Y10ₐ)
      Exo. temp. variation  : $(DstΔTc)
    """

    return jb2008(jd, ϕ_gd, λ, h, F10, F10ₐ, S10, S10ₐ, M10, M10ₐ, Y10, Y10ₐ, DstΔTc)
end

function jb2008(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number,
    F10::Number,
    F10ₐ::Number,
    S10::Number,
    S10ₐ::Number,
    M10::Number,
    M10ₐ::Number,
    Y10::Number,
    Y10ₐ::Number,
    DstΔTc::Number
)
    jd = datetime2julian(instant)
    return jb2008(jd, ϕ_gd, λ, h, F10, F10ₐ, S10, S10ₐ, M10, M10ₐ, Y10, Y10ₐ, DstΔTc)
end

function jb2008(
    jd::Number,
    ϕ_gd::Number,
    λ::Number,
    h::Number,
    F10::Number,
    F10ₐ::Number,
    S10::Number,
    S10ₐ::Number,
    M10::Number,
    M10ₐ::Number,
    Y10::Number,
    Y10ₐ::Number,
    DstΔTc::Number
)
    ########################################################################################
    #                                      Constants                                       #
    ########################################################################################

    T₁    = 183.0       # ............................... Temperature at the lower bound [K]
    z₁    = 90.0        # ................................. Altitude of the lower bound [km]
    zx    = 125.0       # ............................ Altitude of the inflection point [km]
    Rstar = 8314.32     # .. Rstar is the universal gas-constant (mks) [joules / (K . kmol)]
    ρ₁    = 3.46e-6     # ........................................ Density at `z₁` [kg / m³]
    A     = 6.02257e26  # ..................... Avogadro's constant (mks) [molecules / kmol]

    # Assumed sea-level composition.
    Mb₀  = 28.960
    q₀N₂ = 0.78110
    q₀O₂ = 0.20955
    q₀Ar = 9.3400e-3
    q₀He = 1.2890e-5

    # Molecular weights of each specie [kg / kmol].
    MN₂ = 28.0134
    MO₂ = 31.9988
    MO  = 15.9994
    MAr = 39.9480
    MHe =  4.0026
    MH  =  1.00797

    # Thermal diffusion coefficient for each specie.
    α_N₂ =  0.0
    α_O₂ =  0.0
    α_O  =  0.0
    α_Ar =  0.0
    α_He = -0.38
    α_H₂ =  0.0

    # Values used to establish height step sizes in the integration process between 90 km to
    # 105 km, 105 km to 500 km, and above 500km.
    R1 = 0.010
    R2 = 0.025
    R3 = 0.075

    ########################################################################################
    #                                    Preliminaries                                     #
    ########################################################################################

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

    # Compute the local solar time.
    lst = mod((H + π) * 12 / π, 24)

    ########################################################################################
    #                                      Algorithm                                       #
    ########################################################################################

    # == Eq. 2 [1], Eq. 14 [3] =============================================================
    #
    # Nighttime minimum of the global exospheric temperature distribution when the planetary
    # geomagnetic index Kp is zero.
    ΔF10 = F10 - F10ₐ
    ΔS10 = S10 - S10ₐ
    ΔM10 = M10 - M10ₐ
    ΔY10 = Y10 - Y10ₐ

    Wt  = min((F10ₐ / 240)^(1 / 4), 1.0)
    Fsₐ = F10ₐ * Wt + S10ₐ * (1-Wt)
    Tc  = 392.4 + 3.227Fsₐ + 0.298ΔF10 + 2.259ΔS10 + 0.312ΔM10 + 0.178ΔY10

    # == Eq. 15 [3] ========================================================================

    η = abs(ϕ_gd - δs) / 2
    θ = abs(ϕ_gd + δs) / 2

    # == Eq. 16 [3] ========================================================================

    τ = H - 0.64577182 + 0.10471976sin(H + 0.75049158)

    # == Eq. 17 [3] ========================================================================

    m = 2.5
    n = 3.0
    R = 0.31
    C = cos(η)^m
    S = sin(θ)^m

    # NOTE: The original equation in [3] does not have the `abs` as in the source-code of
    # JB2008.
    Tl = Tc * (1 + R * (S + (C - S) * abs(cos(τ / 2))^n))

    # Compute the correction to `Tc` considering the local solar time and latitude.
    ΔTc = _jb2008_ΔTc(F10, lst, ϕ_gd, h)

    # Compute the local exospheric temperature with the geomagnetic storm effect.
    T_exo = Tl + DstΔTc
    T∞    = T_exo + ΔTc

    # == Eq. 9 [3] =========================================================================
    #
    # Temperature at the inflection point `z = 125 km`.
    a  =  444.3807
    b  =    0.02385
    c  = -392.8292
    k  =   -0.0021357
    Tx = a + b * T∞ + c * exp(k * T∞)

    # == Eq. 5 [3] =========================================================================
    #
    # From 90 to 105 km, for a given temperature profile `T[k]`, the density ρ is computed
    # by integrating the barometric equation:
    #
    #                   ┌     ┐
    #                   │  M´ │    M̄´ g
    #   d ln(ρ´) = d ln │ ─── │ - ────── dz ,
    #                   │  T  │    R⋆ T
    #                   └     ┘
    #
    # which can be rewritten as:
    #
    #        ┌       ┐
    #        │  ρ´T  │    M̄´ g
    #   d ln │ ───── │ - ────── dz .
    #        │   M̄´  │    R⋆ T
    #        └       ┘
    #
    # Here, we will compute the following integral:
    #
    #
    #     ╭ z₂
    #     │    M̄´ g
    #     │   ────── dz
    #     │      T
    #     ╯ z₁
    #
    # in which z₁ is the minimum value between `h` and 105 km. The integration will be
    # computed by Newton-Cotes 4th degree method.

    z₂ = min(h, 105.0)

    int, z₂ = _jb2008_∫(z₁, z₂, R1, Tx, T∞, _jb2008_δf1)

    Mb₁ = _jb2008_mean_molecular_mass(z₁)
    Tl₁ = _jb2008_temperature(z₁, Tx, T∞)
    Mb₂ = _jb2008_mean_molecular_mass(z₂)
    Tl₂ = _jb2008_temperature(z₂, Tx, T∞)

    # `Mbj` and `Tlj` contains, respectively, the mean molecular mass and local temperature
    # at the boundary of the integration interval.
    #
    # The factor 1000 converts `Rstar` to the appropriate units.
    ρ = ρ₁ * (Mb₂ / Mb₁) * (Tl₁ / Tl₂) * exp(-1000 * int / Rstar)

    # == Eq. 2 [3] =========================================================================

    NM = A * ρ
    N  = NM / Mb₂

    # == Eq. 3 [3] =========================================================================

    log_nN₂ = log(q₀N₂ * NM / Mb₀)
    log_nAr = log(q₀Ar * NM / Mb₀)
    log_nHe = log(q₀He * NM / Mb₀)

    # == Eq. 4 [3] =========================================================================

    log_nO₂ = log(NM / Mb₀ * (1 + q₀O₂) - N)
    log_nO  = log(2 * (N - NM / Mb₀))
    log_nH  = 0.0

    Tz = 0.0

    if h <= 105
        Tz = Tl₂
        log_nH₂ = log_nHe - 25

    else
        # == Eq.6 [3] ======================================================================
        #
        # From 100 km to 500 km, neglecting the thermal diffusion coefficient,
        # the eq. 6 in [3] can be written as:
        #
        #        ┌            ┐
        #        │      1 + αᵢ│      M̄´ g
        #   d ln │n . T       │ = - ────── dz ,
        #        │            │      R* T
        #        └            ┘
        #
        # where `n` is number density for the i-th specie. This equations must be integrated
        # for each specie we are considering.
        #
        # Here, we will compute the following integral:
        #
        #   ╭ z₃
        #   │    g
        #   │   ─── dz
        #   │    T
        #   ╯ z₂
        #
        # in which z₃ is the minimum value between `h` and 500 km. The integration will be
        # computed by Newton-Cotes 4th degree method.

        z₃ = min(h, 500.0)

        int₁, z₃ = _jb2008_∫(z₂, z₃, R1, Tx, T∞, _jb2008_δf2)

        Tl₃ = _jb2008_temperature(z₃, Tx, T∞)

        # If `h` is lower than 500 km, then keep integrating until 500 km due to the
        # hydrogen. Otherwise, continue the integration until `h`. Hence, we will compute
        # the following integral:
        #
        #   ╭ z₄
        #   │    g
        #   │   ─── dz
        #   │    T
        #   ╯ z₃
        #
        # in which z₄ is the maximum value between `h` and 500 km. The integration will be
        # computed by Newton-Cotes 4th degree method.

        z₄ = max(h, 500.0)

        int₂, z₄ = _jb2008_∫(z₃, z₄, (h <= 500) ? R2 : R3, Tx, T∞, _jb2008_δf2)

        Tl₄ = _jb2008_temperature(z₄, Tx, T∞)

        if h <= 500
            Tz        = Tl₃
            log_TfoTi = log(Tl₃ / Tl₂)
            H_sign    = +1

            #           g
            # goRT =  ─────
            #         R⋆  T
            goRT = 1000 * int₁ / Rstar
        else
            Tz        = Tl₄
            log_TfoTi = log(Tl₄ / Tl₂)
            H_sign    = -1

            #           g
            # goRT =  ─────
            #         R⋆  T
            goRT = 1000 * (int₁ + int₂) / Rstar
        end

        log_nN₂ += -(1 + α_N₂) * log_TfoTi - goRT * MN₂
        log_nO₂ += -(1 + α_O₂) * log_TfoTi - goRT * MO₂
        log_nO  += -(1 + α_O ) * log_TfoTi - goRT * MO
        log_nAr += -(1 + α_Ar) * log_TfoTi - goRT * MAr
        log_nHe += -(1 + α_He) * log_TfoTi - goRT * MHe

        # == Eq. 7 [3] =====================================================================
        #
        # The equation of [3] will be converted from `log₁₀` to `log`.  Furthermore, the
        # units will be converted to [1 / cm³] to [1 / m³].

        log₁₀_nH_500km = @evalpoly(log10(T∞), 73.13, -39.40, 5.5)
        log_nH_500km   = log(10) * (log₁₀_nH_500km + 6)
        #                                            ^
        #                                            |
        #                       This factor converts from [1 / cm³] to [1 / m³].

        # Compute the hydrogen density based on the value at 500km.
        log_nH = log_nH_500km + H_sign * (log(Tl₄ / Tl₃) + 1000 / Rstar * int₂ * MH)
    end

    # == Eq. 24 [3] - Seasonal-latitudinal variation =======================================
    #
    # TODO: The term related to the year in the source-code of JB 2008 is different from
    # [3]. This must be verified.

    #       | Modified jd  |
    Φ = mod((jd - 2400000.5 - 36204 ) / 365.2422, 1)

    Δlog₁₀ρ = 0.02 * (h - 90) * sign(ϕ_gd) * exp(-0.045(h - 90)) * sin(ϕ_gd)^2 * sin(2π * Φ + 1.72 )

    # Convert from `log10` to `log`.
    Δlogρ = log(10) * Δlog₁₀ρ

    # == Eq. 23 [3] - Semiannual variation =================================================

    if h < 2000
        # Compute the year given the selected Julian Day.
        year, month, day, = jd_to_date(jd)

        # Compute the day of the year.
        doy = jd - date_to_jd(year, 1, 1, 0, 0, 0) + 1

        # Use the new semiannual model from [1].
        Fz, Gz, Δsalog₁₀ρ = _jb2008_semiannual(doy, h, F10ₐ, S10ₐ, M10ₐ)

        (Fz < 0) && (Δsalog₁₀ρ = 0.0)

        # Convert from `log10` to `log`.
        Δsalogρ = log(10) * Δsalog₁₀ρ

    else
        Δsalogρ = 0.0
    end

    # Compute the total variation.
    Δρ = Δlogρ + Δsalogρ

    # Apply to the number densities.
    log_nN₂ += Δρ
    log_nO₂ += Δρ
    log_nO  += Δρ
    log_nAr += Δρ
    log_nHe += Δρ
    log_nH  += Δρ

    # Compute the high altitude exospheric density correction factor.
    #
    # The source-code computes this correction after computing the mass density and mean
    # molecular weight. Hence, the correction is applied only to the total density. However,
    # we will apply the correction to each specie.
    #
    # TODO: Verify if this is reasonable.

    log_FρH  = log(_jb2008_high_altitude(h, F10ₐ))
    log_nN₂ += log_FρH
    log_nO₂ += log_FρH
    log_nO  += log_FρH
    log_nAr += log_FρH
    log_nHe += log_FρH
    log_nH  += log_FρH

    # Compute the mass density and mean molecular weight and convert number density logs
    # from natural to common.
    sum_n  = 0.0
    sum_mn = 0.0

    for (log_n, M) in (
        (log_nN₂, MN₂),
        (log_nO₂, MO₂),
        (log_nO, MO),
        (log_nAr, MAr),
        (log_nHe, MHe),
        (log_nH, MH)
    )
        n = exp(log_n)
        sum_n  += n
        sum_mn += n * M
    end

    nN₂ = exp(log_nN₂)
    nO₂ = exp(log_nO₂)
    nO  = exp(log_nO)
    nAr = exp(log_nAr)
    nHe = exp(log_nHe)
    nH  = exp(log_nH)
    ρ   = sum_mn / A

    # Create and return the output structure.
    return JB2008Output{Float64}(
        ρ,
        Tz,
        T_exo,
        nN₂,
        nO₂,
        nO,
        nAr,
        nHe,
        nH,
    )
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

#   _jb2008_gravity(z::Number) -> Float64
#
# Compute the gravity [m / s] at altitude `z` [km] according to the model Jacchia 1971 [3].
function _jb2008_gravity(z::Number)
    # Mean Earth radius [km].
    Re = 6356.776

    # Gravity at Earth surface [m/s²].
    g₀ = 9.80665

    # Gravity at desired altitude [m/s²].
    return g₀ / (1 + z / Re)^2
end

#   _jb2008_high_altitude(h::Number, F10ₐ::Number) -> Float64
#
# Compute the high altitude exospheric density correction factor in altitude `h` [km] and
# the averaged 10.7-cm solar flux `F10ₐ` [sfu], which must be obtained using 81-day centered
# on input time.
#
# This function uses the model in Section 6.2 of [2].
function _jb2008_high_altitude(h::Number, F10ₐ::Number)
    # Auxiliary variables.
    C = _JB2008_CHT

    # Compute the high-altitude density correction.
    FρH = 1.0

    @inbounds if 1000 <= h <= 1500
        z = (h - 1000) / 500

        # In this case, the `F₁₅₀₀` is the density factor at 1500 km.
        F₁₅₀₀ = C[1] + C[2] * F10ₐ + 1500 * (C[3] + C[4] * F10ₐ)

        ∂F₁₅₀₀_∂z = 500 * (C[3] + C[4] * F10ₐ)

        # In [2], there is an error when computing this coefficients (eq. 11). The partial
        # derivative has the `500` factor and the coefficients have other `500` factor that
        # **is not** present in the source-code.
        c0  = +1
        c1  = +0
        c2  = +3F₁₅₀₀ - ∂F₁₅₀₀_∂z - 3
        c3  = -2F₁₅₀₀ + ∂F₁₅₀₀_∂z + 2
        FρH = @evalpoly(z, c0, c1, c2, c3)

    elseif h > 1500
        FρH = C[1] + C[2] * F10ₐ + C[3] * h + C[4] * h * F10ₐ

    end

    return FρH
end

#   _jb2008_M(z::R) where R -> Float64
#
# Compute the mean molecular mass at altitude `z` [km] using the empirical profile in eq. 1
# [3].
function _jb2008_mean_molecular_mass(z::Number)
    !(90 <= z < 105.1) &&
        @warn "The empirical model for the mean molecular mass is valid only for 90 <= z <= 105 km."

    M = @evalpoly(
        z - 100,
        +28.15204,
         -0.085586,
         +1.2840e-4,
         -1.0056e-5,
         -1.0210e-5,
         +1.5044e-6,
         +9.9826e-8
    )

    return M
end

#   _jb2008_temperature(z::Number, Tx::Number, T∞::Number)
#
# Compute the temperature [K] at height `z` [km] given the temperature `Tx` [K] at the
# inflection point, and the exospheric temperature `T∞` [K] according to the theory of the
# model Jacchia 1971 [3].
#
# The inflection point is considered to by `z = 125 km`.
function _jb2008_temperature(z::Number, Tx::Number, T∞::Number)
    # == Constants =========================================================================

    T₁  = 183      # .................................... Temperature at the lower bound [K]
    z₁  = 90       # ...................................... Altitude of the lower bound [km]
    zx  = 125      # ................................. Altitude of the inflection point [km]
    Δz₁ = z₁ - zx

    # == Check the Parameters ==============================================================

    (z  < z₁) && throw(ArgumentError("The altitude must not be lower than $(z₁) km."))
    (T∞ < 0)  && throw(ArgumentError("The exospheric temperature must be positive."))

    # Compute the temperature gradient at the inflection point.
    Gx = 1.9 * (Tx - T₁) / (zx - z₁)

    # == Compute the Temperature at the Desire Altitude ====================================

    Δz = z - zx

    if z <= zx
        c₁ = Gx
        c₂ = 0
        c₃ =-5.1 / (5.7 * Δz₁^2) * Gx
        c₄ = 0.8 / (1.9 * Δz₁^3) * Gx
        T  = @evalpoly(Δz, Tx, c₁, c₂, c₃, c₄)
    else
        A  = 2 * (T∞ - Tx) / π
        T  = Tx + A * atan(Gx / A * Δz * (1 + 4.5e-6 * Δz^2.5))
    end

    return T
end


#   _jb2008_δf1(z::Number, Tx::Number, T∞::Number) -> Float64
#
# Auxiliary function to compute the integrand in `_jb2008_∫`.
function _jb2008_δf1(z::Number, Tx::Number, T∞::Number)
    Mb = _jb2008_mean_molecular_mass(z)
    Tl = _jb2008_temperature(z, Tx, T∞)
    g  = _jb2008_gravity(z)

    return Mb * g / Tl
end

#   _jb2008_δf2(z, Tx, T∞)
#
# Auxiliary function to compute the integrand in `_jb2008_∫`.
function _jb2008_δf2(z::Number, Tx::Number, T∞::Number)
    Tl = _jb2008_temperature(z, Tx, T∞)
    g  = _jb2008_gravity(z)

    return g / Tl
end

#   _jb2008_∫(z₀::Number, z₁::Number, R::Number, Tx::Number, T∞::Number, δf::Function) -> Float64, Float64
#
# Compute the integral of the function `δf` between `z₀` and `z₁` using the Newton-Cotes 4th
# degree method. `R` is a number that defines the step size, `Tx` is the temperature at the
# inflection point, and `T∞` is the exospheric temperature.
#
# The integrand function `δf` must be `_jb2008_δf1` or `_jb2008_δf2`.
#
# # Returns
#
# - `Float64`: Value of the integral.
# - `Float64`: Last value of `z` used in the numerical algorithm.
function _jb2008_∫(
    z₀::Number,
    z₁::Number,
    R::Number,
    Tx::Number,
    T∞::Number,
    δf::Function
)
    # Compute the number of integration steps.
    #
    # This is computed so that `z₂ = z₁*(zr)^n`. Hence, `zr` is the factor that defines the
    # size of each integration interval.
    al = log(z₁ / z₀)
    n  = floor(Int64, al / R) + 1
    zr = exp(al / n)

    # Initialize the integration auxiliary variables.
    zi₁ = z₀
    zj  = 0.0
    Mbj = 0.0
    Tlj = 0.0

    # Variable to store the integral from `z₀` to `z₁`.
    int = 0.0

    # For each integration step, use the Newton-Cotes 4th degree formula to integrate
    # (Boole's rule).
    @inbounds for i in 1:n
        zi₀   = zi₁             # ............... The beginning of the i-th integration step
        zi₁   = zr * zi₁        # ..................... The end of the i-th integration step
        Δz    = (zi₁ - zi₀) / 4 # ....................... Step for the i-th integration step

        # Compute the Newton-Cotes 4th degree sum.
        zj     = zi₀
        int_i  = 14 // 45 * δf(zj, Tx, T∞)

        zj    += Δz
        int_i += 64 // 45 * δf(zj, Tx, T∞)

        zj    += Δz
        int_i += 24 // 45 * δf(zj, Tx, T∞)

        zj    += Δz
        int_i += 64 // 45 * δf(zj, Tx, T∞)

        zj    += Δz
        int_i += 14 // 45 * δf(zj, Tx, T∞)

        # Accumulate the sum.
        int   += int_i * Δz
    end

    return int, zj
end

#
#   _jb2008_semiannual(doy::Number, h::Number, F10ₐ::Number, S10ₐ::Number, M10ₐ::Number)
#
# Compute the semiannual density variation considering the JB2008 model [1].
#
# # Arguments
#
# - `doy::Number`: Day of the year + fraction of the day.
# - `h::Number`: Height [km].
# - `F10ₐ::Number`: Averaged 10.7-cm flux, which must use a 81-day window centered on
#   the input time [10⁻²² W / (M² Hz)].
# - `S10ₐ`: EUV 81-day averaged centered index.
# - `M10ₐ`: MG2 81-day averaged centered index.
#
# # Returns
#
# - `Float64`: Semiannual F(z) height function.
# - `Float64`: Semiannual G(t) yearly periodic function.
# - `Float64`: Semiannual variation of the density `Δsalog₁₀ρ`.
function _jb2008_semiannual(
    doy::Number,
    h::Number,
    F10ₐ::Number,
    S10ₐ::Number,
    M10ₐ::Number
)
    # Auxiliary variables.
    B = _JB2008_FZM
    C = _JB2008_GTM
    z = h / 1000
    ω = 2π * (doy - 1) / 365  # .............................................. See eq. 4 [2]

    # Compute the new 81-day centered solar index for F(z) according to eq. 4 [1].
    Fsmjₐ = 1F10ₐ - 0.7S10ₐ - 0.04M10ₐ

    # Compute the semiannual F(z) height function according to eq. 5 [1].
    Fz = B[1] + B[2] * Fsmjₐ + B[3] * z * Fsmjₐ + B[4] * z^2 * Fsmjₐ + B[5] * z * Fsmjₐ^2

    # Compute the new 81-day centered solar index for G(t) according to eq. 6 [1].
    Fsmₐ = 1F10ₐ - 0.75S10ₐ - 0.37M10ₐ

    # Compute the semiannual G(t) yearly periodic function according to eq. 7
    # [1].
    sω,  cω  = sincos(1ω)
    s2ω, c2ω = sincos(2ω)

    Gt = (C[1] + C[2] * sω + C[3] * cω + C[4] * s2ω + C[5]  * c2ω) +
         (C[6] + C[7] * sω + C[8] * cω + C[9] * s2ω + C[10] * c2ω) * Fsmₐ

    Fz = max(Fz, 1e-6)

    Δsalog₁₀ρ = Fz * Gt

    return Fz, Gt, Δsalog₁₀ρ
end

#   _jb2008_ΔTc(F10::Number, lst::Number, ϕ_gd::Number, h::Number)
#
# Compute the correction in the `Tc` for Jacchia-Bowman model.
#
# This correction is mention in [2]. However, the equations do not seem to match those in
# the source-code. The ones implemented here are exactly the same as in the source-code.
#
# # Arguments
#
# - `F10::Number`: F10.7 flux.
# - `lst::Number`: Local solar time (0 - 24) [hr].
# - `ϕ_gd::Number`: Geodetic latitude [rad].
# - `h::Number`: Altitude [km].
#
# # Returns
#
# - `Float64`: The correction `ΔTc` [K].
function _jb2008_ΔTc(F10::Number, lst::Number, ϕ_gd::Number, h::Number)
    # Auxiliary variables according to [2, p.  784].
    B  = _JB2008_B
    C  = _JB2008_C
    F  = (F10 - 100) / 100
    θ  = lst / 24
    θ² = θ^2
    θ³ = θ^3
    θ⁴ = θ^4
    θ⁵ = θ^5
    cϕ = cos(ϕ_gd)

    ΔTc = 0.0

    # Compute the temperature variation given the altitude.
    @inbounds if 120 <= h <= 200
        ΔTc200 = C[17] +
                 C[18] * θ  * cϕ +
                 C[19] * θ² * cϕ +
                 C[20] * θ³ * cϕ +
                 C[21] * F  * cϕ +
                 C[22] * θ  * F * cϕ +
                 C[23] * θ² * F * cϕ

        ΔTc200Δz = C[ 1] +
                   B[ 2] * F +
                   C[ 3] * θ  * F +
                   C[ 4] * θ² * F +
                   C[ 5] * θ³ * F +
                   C[ 6] * θ⁴ * F +
                   C[ 7] * θ⁵ * F +
                   C[ 8] * θ  * cϕ +
                   C[ 9] * θ² * cϕ +
                   C[10] * θ³ * cϕ +
                   C[11] * θ⁴ * cϕ +
                   C[12] * θ⁵ * cϕ +
                   C[13] * cϕ +
                   C[14] * F  * cϕ +
                   C[15] * θ  * F * cϕ +
                   C[16] * θ² * F * cϕ

        zp  = (h - 120) / 80
        ΔTc = (3ΔTc200 - ΔTc200Δz) * zp^2 + (ΔTc200Δz - 2ΔTc200) * zp^3

    elseif 200 < h <= 240
        H = (h - 200) / 50

        ΔTc = C[ 1] * H +
              B[ 2] * F  * H +
              C[ 3] * θ  * F  * H +
              C[ 4] * θ² * F  * H +
              C[ 5] * θ³ * F  * H +
              C[ 6] * θ⁴ * F  * H +
              C[ 7] * θ⁵ * F  * H +
              C[ 8] * θ  * cϕ * H +
              C[ 9] * θ² * cϕ * H +
              C[10] * θ³ * cϕ * H +
              C[11] * θ⁴ * cϕ * H +
              C[12] * θ⁵ * cϕ * H +
              C[13] * cϕ * H +
              C[14] * F  * cϕ * H +
              C[15] * θ  * F  * cϕ * H +
              C[16] * θ² * F  * cϕ * H +
              C[17] +
              C[18] * θ  * cϕ +
              C[19] * θ² * cϕ +
              C[20] * θ³ * cϕ +
              C[21] * F  * cϕ +
              C[22] * θ  * F  * cϕ +
              C[23] * θ² * F  * cϕ

    elseif 240 < h <= 300
        H = 40 / 50

        aux1 = C[ 1] * H +
               B[ 2] * F  * H +
               C[ 3] * θ  * F  * H +
               C[ 4] * θ² * F  * H +
               C[ 5] * θ³ * F  * H +
               C[ 6] * θ⁴ * F  * H +
               C[ 7] * θ⁵ * F  * H +
               C[ 8] * θ  * cϕ * H +
               C[ 9] * θ² * cϕ * H +
               C[10] * θ³ * cϕ * H +
               C[11] * θ⁴ * cϕ * H +
               C[12] * θ⁵ * cϕ * H +
               C[13] * cϕ * H +
               C[14] * F  * cϕ * H +
               C[15] * θ  * F  * cϕ * H +
               C[16] * θ² * F  * cϕ * H +
               C[17] +
               C[18] * θ  * cϕ +
               C[19] * θ² * cϕ +
               C[20] * θ³ * cϕ +
               C[21] * F  * cϕ +
               C[22] * θ  * F  * cϕ +
               C[23] * θ² * F  * cϕ

        aux2 = C[ 1] +
               B[ 2] * F +
               C[ 3] * θ  * F +
               C[ 4] * θ² * F +
               C[ 5] * θ³ * F +
               C[ 6] * θ⁴ * F +
               C[ 7] * θ⁵ * F +
               C[ 8] * θ  * cϕ +
               C[ 9] * θ² * cϕ +
               C[10] * θ³ * cϕ +
               C[11] * θ⁴ * cϕ +
               C[12] * θ⁵ * cϕ +
               C[13] * cϕ +
               C[14] * F  * cϕ +
               C[15] * θ  * F * cϕ +
               C[16] * θ² * F * cϕ

        H = 300 / 100

        ΔTc300 = B[ 1] +
                 B[ 2] * F +
                 B[ 3] * θ  * F +
                 B[ 4] * θ² * F +
                 B[ 5] * θ³ * F +
                 B[ 6] * θ⁴ * F +
                 B[ 7] * θ⁵ * F +
                 B[ 8] * θ  * cϕ +
                 B[ 9] * θ² * cϕ +
                 B[10] * θ³ * cϕ +
                 B[11] * θ⁴ * cϕ +
                 B[12] * θ⁵ * cϕ +
                 B[13] * H  * cϕ +
                 B[14] * θ  * H * cϕ +
                 B[15] * θ² * H * cϕ +
                 B[16] * θ³ * H * cϕ +
                 B[17] * θ⁴ * H * cϕ +
                 B[18] * θ⁵ * H * cϕ +
                 B[19] * cϕ

        ΔTc300Δz = B[13] * cϕ +
                   B[14] * θ  * cϕ +
                   B[15] * θ² * cϕ +
                   B[16] * θ³ * cϕ +
                   B[17] * θ⁴ * cϕ +
                   B[18] * θ⁵ * cϕ

        aux3 = 3ΔTc300 - ΔTc300Δz - 3aux1 - 2aux2
        aux4 = ΔTc300 - aux1 - aux2 - aux3
        zp   = (h - 240) / 60
        ΔTc  = @evalpoly(zp, aux1, aux2, aux3, aux4)

    elseif 300 < h <= 600
         H  = h / 100

        ΔTc = B[ 1] +
              B[ 2] * F +
              B[ 3] * θ  * F +
              B[ 4] * θ² * F +
              B[ 5] * θ³ * F +
              B[ 6] * θ⁴ * F +
              B[ 7] * θ⁵ * F +
              B[ 8] * θ  * cϕ +
              B[ 9] * θ² * cϕ +
              B[10] * θ³ * cϕ +
              B[11] * θ⁴ * cϕ +
              B[12] * θ⁵ * cϕ +
              B[13] * H  * cϕ +
              B[14] * θ  * H * cϕ +
              B[15] * θ² * H * cϕ +
              B[16] * θ³ * H * cϕ +
              B[17] * θ⁴ * H * cϕ +
              B[18] * θ⁵ * H * cϕ +
              B[19] * cϕ

    elseif 600 < h <= 800
        zp   = (h - 600) / 100
        hp   = 600 / 100

        aux1 = B[ 1] +
               B[ 2] * F +
               B[ 3] * θ  * F +
               B[ 4] * θ² * F +
               B[ 5] * θ³ * F +
               B[ 6] * θ⁴ * F +
               B[ 7] * θ⁵ * F +
               B[ 8] * θ  * cϕ +
               B[ 9] * θ² * cϕ +
               B[10] * θ³ * cϕ +
               B[11] * θ⁴ * cϕ +
               B[12] * θ⁵ * cϕ +
               B[13] * hp * cϕ +
               B[14] * θ  * hp * cϕ +
               B[15] * θ² * hp * cϕ +
               B[16] * θ³ * hp * cϕ +
               B[17] * θ⁴ * hp * cϕ +
               B[18] * θ⁵ * hp * cϕ +
               B[19] * cϕ

        aux2 = B[13] * cϕ +
               B[14] * θ  * cϕ +
               B[15] * θ² * cϕ +
               B[16] * θ³ * cϕ +
               B[17] * θ⁴ * cϕ +
               B[18] * θ⁵ * cϕ

        aux3 = -(3aux1 + 4aux2) / 4
        aux4 =  ( aux1 +  aux2) / 4
        ΔTc  = @evalpoly(zp, aux1, aux2, aux3, aux4)
    end

    return ΔTc
end
