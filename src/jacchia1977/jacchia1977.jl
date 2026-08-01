## Description #############################################################################
#
# The Jacchia 1977 Atmospheric Model.
#
## References ##############################################################################
#
# [1] Jacchia, L. G (1977). Thermospheric temperature, density and composition: New models.
#     SAO Special Report #375.
#
# [2] de Matos, B. S., Carrara, V (1985-1987). Fortran implementation of the Jacchia 1977
#     model (routines IDYMOS, ISDAMO, ISMADE, DIVARI, GEOACI, IMOWEI, TEMLO, SEALAT, and
#     SEMIAN). INPE, São José dos Campos, BR. Available in:
#
#         https://github.com/jacobwilliams/INPE-atmosphere-models
#
############################################################################################

export jacchia1977

"""
    jacchia1977(
        instant::DateTime,
        ϕ_gd::Number,
        λ::Number,
        h::Number[, F10::Number, F10ₐ::Number, Kp::Number];
        kwargs...
    ) -> Jacchia1977Output
    jacchia1977(
        jd::Number,
        ϕ_gd::Number,
        λ::Number,
        h::Number[, F10::Number, F10ₐ::Number, Kp::Number];
        kwargs...
    ) -> Jacchia1977Output

Compute the atmospheric density using the Jacchia 1977 model.

Unlike the Jacchia-Roberts 1971 model, the Jacchia 1977 model does not have a closed-form
solution. Hence, the barometric and diffusion equations are numerically integrated here,
making this model considerably slower than [`jr1971`](@ref).

If we omit all space indices, the system tries to obtain them automatically for the
selected day `jd` or `instant` using the prescriptions in the report [1]: the daily flux is
evaluated with a lag that depends on the solar hour angle, the averaged flux is a
Gaussian-weighted mean with a standard width of 71 days centered on the input time, and the
Kp is delayed by an interval that depends on the geomagnetic latitude. However, the indices
must be already initialized using the function `SpaceIndices.init()`. Notice that the
Gaussian-weighted mean requires the indices to be available up to 213 days after the input
time.

The function throws an `ArgumentError` if the altitude `h` is outside the interval
[90, 2000] km.

# Arguments

- `jd::Number`: Julian day to compute the model.
- `instant::DateTime`: Instant to compute the model represented using `DateTime`.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `h::Number`: Altitude [m].
- `F10::Number`: 10.7-cm solar flux [sfu], evaluated with the lag prescribed in [1] (from
    0.9 day to 1.6 days, depending on the local solar time).
- `F10ₐ::Number`: 10.7-cm averaged solar flux, Gaussian-weighted mean centered on the input
    time with a standard width of 71 days [sfu].
- `Kp::Number`: Kp geomagnetic index, delayed by 0.1 day to 0.3 day depending on the
    geomagnetic latitude, as prescribed in [1].

# Keywords

- `verbose::Val`: Set to `Val(true)` to emit debug messages related to the automatic space
    index fetching, or to `Val(false)` to suppress them. Notice that this keyword must be a
    `Val` object, not a `Bool`, and it is only available in the methods that fetch the
    space indices automatically.
    (**Default**: `Val(true)`)

# Returns

- `Jacchia1977Output`: Structure containing the results obtained from the model. Its
    element type is the promotion of the types of the numeric inputs.

# References

- **[1]** Jacchia, L. G (1977). *Thermospheric temperature, density and composition: New
    models*. SAO Special Report #375.
- **[2]** de Matos, B. S., Carrara, V (1985-1987). *Fortran implementation of the Jacchia
    1977 model*. INPE, São José dos Campos, BR.
"""
function jacchia1977(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number;
    verbose::Val{verbosity} = Val(true),
) where {verbosity}
    return jacchia1977(datetime2julian(instant), ϕ_gd, λ, h; verbose = verbose)
end

function jacchia1977(
    jd::Number, ϕ_gd::Number, λ::Number, h::Number; verbose::Val{verbosity} = Val(true)
) where {verbosity}
    # == Daily F10.7 With the Solar Hour Angle Dependent Lag, Eq. 23 [1] ===================

    # Compute the Sun position represented in the inertial reference frame (MOD).
    s_i = sun_position_mod(jd)

    # Compute the Sun right ascension [rad].
    Ωs = atan(s_i[2], s_i[1])

    # Compute the right ascension of the selected location w.r.t. the inertial reference
    # frame.
    Ωp = λ + jd_to_gmst(jd)

    # Hour angle of the Sun at the selected location [rad].
    H = Ωp - Ωs

    # Lag of the daily flux [days].
    Δt = 1.26 + 0.37 * sin(H - deg2rad(92))

    F10 = space_index(Val(:F10obs), jd - Δt)

    # == Gaussian-Weighted Averaged F10.7, Eqs. 21-22 [1] ==================================

    # The report recommends a Gaussian weight with a standard width of three solar
    # rotations (71 days). We truncate the window at three standard widths.
    Σw   = 0.0
    Σw_F = 0.0

    for k in -213:213
        w    = exp(-(k / 71)^2)
        Σw   += w
        Σw_F += w * space_index(Val(:F10obs), jd + k)
    end

    F10ₐ = Σw_F / Σw

    # == Kp Delayed by the Geomagnetic Latitude Dependent Lag, Eq. 30 [1] ==================

    # Compute the sine of the geomagnetic (invariant) latitude using the dipole
    # approximation in [2].
    sin_ϕᵢ = 0.9792 * sin(ϕ_gd) + 0.2028 * cos(ϕ_gd) * cos(λ - deg2rad(291))

    # Lag of the geomagnetic index [days].
    τ = 0.1 + 0.2 * (1 - sin_ϕᵢ^2)

    instant_delayed = julian2datetime(jd - τ)
    Kp_vect = space_index(Val(:Kp), instant_delayed)

    # Get the number of seconds elapsed since the beginning of the day of the delayed
    # instant and select the related 3-hour interval.
    day = Date(instant_delayed) |> DateTime
    Δts = Dates.value(instant_delayed - day) / 1000
    id  = clamp(floor(Int, Δts / 10_800) + 1, 1, 8)
    Kp  = Kp_vect[id]

    verbosity && @debug """
    Jacchia 1977 - Fetched Space Indices
      Lagged daily F10.7      : $(F10) sfu (lag = $(Δt) days)
      Gaussian averaged F10.7 : $(F10ₐ) sfu
      Delayed Kp              : $(Kp) (lag = $(τ) days)
    """

    return jacchia1977(jd, ϕ_gd, λ, h, F10, F10ₐ, Kp)
end

function jacchia1977(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number,
    F10::Number,
    F10ₐ::Number,
    Kp::Number,
)
    return jacchia1977(datetime2julian(instant), ϕ_gd, λ, h, F10, F10ₐ, Kp)
end

function jacchia1977(
    jd::JT, ϕ_gd::PT, λ::LT, h::HT, F10::FT, F10ₐ::FT2, Kp::KT
) where {
    JT <: Number,
    PT <: Number,
    LT <: Number,
    HT <: Number,
    FT <: Number,
    FT2 <: Number,
    KT <: Number,
}
    RT = float(promote_type(JT, PT, LT, HT, FT, FT2, KT))

    # == Preliminaries =====================================================================

    # Convert the altitude from [m] to [km] and check the bounds.
    z = h / 1000

    !(90 <= z <= 2000) &&
        throw(ArgumentError("The altitude must be between 90 km and 2000 km."))

    # Compute the Sun position represented in the inertial reference frame (MOD).
    s_i = sun_position_mod(jd)

    # Compute the Sun declination [rad].
    δs = atan(s_i[3], √(s_i[1] * s_i[1] + s_i[2] * s_i[2]))

    # Compute the Sun right ascension [rad].
    Ωs = atan(s_i[2], s_i[1])

    # Compute the right ascension of the selected location w.r.t. the inertial reference
    # frame.
    Ωp = λ + jd_to_gmst(jd)

    # Fraction of the tropic year starting on January 1st, as in [2] (the epoch is the
    # modified Julian date referred to 1950.0).
    Φ = mod((jd - 2433282.5) / 365.2422, 1)

    return _jacchia1977_dynamic(RT(z), ϕ_gd, Ωp, Ωs, δs, λ, Φ, F10, F10ₐ, Kp)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _jacchia1977_dynamic(
        z::Number,
        ϕ::Number,
        Ωp::Number,
        Ωs::Number,
        δs::Number,
        λ::Number,
        Φ::Number,
        F10::Number,
        F10ₐ::Number,
        Kp::Number
    ) -> Jacchia1977Output

Compute the Jacchia 1977 dynamic model (routine ISDAMO of [2]).

# Arguments

- `z::Number`: Altitude [km].
- `ϕ::Number`: Latitude [rad].
- `Ωp::Number`: Right ascension of the selected location [rad].
- `Ωs::Number`: Right ascension of the Sun [rad].
- `δs::Number`: Declination of the Sun [rad].
- `λ::Number`: Longitude of the selected location [rad].
- `Φ::Number`: Fraction of the tropic year starting on January 1st [-].
- `F10::Number`: 10.7-cm solar flux [sfu].
- `F10ₐ::Number`: 10.7-cm averaged solar flux [sfu].
- `Kp::Number`: Kp geomagnetic index [-].

# Returns

- `Jacchia1977Output`: Structure containing the results obtained from the model. Its
    element type is the promotion of the types of the numeric inputs.
"""
function _jacchia1977_dynamic(
    z::Number,
    ϕ::Number,
    Ωp::Number,
    Ωs::Number,
    δs::Number,
    λ::Number,
    Φ::Number,
    F10::Number,
    F10ₐ::Number,
    Kp::Number,
)
    Mi = _JACCHIA1977_CONSTANTS.Mi
    Av = _JACCHIA1977_CONSTANTS.Av

    RT = float(
        promote_type(
            typeof(z),
            typeof(ϕ),
            typeof(Ωp),
            typeof(Ωs),
            typeof(δs),
            typeof(λ),
            typeof(Φ),
            typeof(F10),
            typeof(F10ₐ),
            typeof(Kp),
        ),
    )

    # == Static Model at the Mean Exospheric Temperature, Eq. 20 [1] =======================

    T½ = 5.48 * F10ₐ^RT(0.8) + 101.8 * F10^RT(0.4)

    _, M̄, _ = _jacchia1977_static(T½, z)

    # Local temperature from the static temperature profile.
    Tz = _jacchia1977_temperature(z, _jacchia1977_profile_params(T½))

    # == Diurnal Variation, Eqs. 24-27 [1] =================================================

    ad, ac, Θ_H = _jacchia1977_diurnal(T½, Ωp, Ωs, δs, ϕ, z, M̄)

    # == Geomagnetic Variation, Eqs. 28-35 [1] =============================================

    dn = _jacchia1977_geomagnetic(Θ_H, Kp, ϕ, λ, z)
    ad = ad .+ dn .- ac

    # == Seasonal-Latitudinal Variation, Eqs. 36-39 [1] ====================================

    ad = ad .+ _jacchia1977_seasonal_latitudinal(Φ, δs, ϕ, z)

    # == Semiannual Variation, Eqs. 40-44 [1] ==============================================

    ad = ad .+ _jacchia1977_semiannual(Φ, z)

    # == Assemble the Output ===============================================================

    n = 10 .^ ad
    ρ = sum(n .* Mi) / Av

    return Jacchia1977Output{RT}(ρ, Tz, T½, n[3], n[2], n[5], n[4], n[1], n[6])
end

"""
    _jacchia1977_profile_params(T∞::Number) -> NTuple{7, Number}

Compute the parameters of the temperature profile (eqs. 1 to 4 of [1]) related to the
exospheric temperature `T∞` [K].

See also: [`_jacchia1977_temperature`](@ref)
"""
function _jacchia1977_profile_params(T∞::Number)
    T₀ = _JACCHIA1977_CONSTANTS.T₀
    z₀ = _JACCHIA1977_CONSTANTS.z₀
    zx = _JACCHIA1977_CONSTANTS.zx

    # Temperature at the inflection point, eq. 1 [1].
    aux = 0.0045 * (T∞ - T₀)
    Tx  = T₀ + 110.5 * atanh(aux / √(1 + aux * aux))

    # Temperature gradient at the inflection point, eq. 2 [1].
    Gx = 1.9 * (Tx - T₀) / (zx - z₀)

    c₁ = 2 * (Tx - T₀) / π
    c₂ = Gx / c₁
    c₃ = 1.7 * c₂
    c₄ = 2 * (T∞ - Tx) / π
    c₅ = Gx / c₄
    c₆ = 5.5e-5 * c₅

    return (c₁, c₂, c₃, c₄, c₅, c₆, Tx)
end

"""
    _jacchia1977_temperature(z::Number, c::NTuple{7, T}) where {T<:Number} -> Number

Compute the temperature [K] at the altitude `z` [km] using the profile parameters `c`
obtained from [`_jacchia1977_profile_params`](@ref) (eqs. 3 and 4 of [1]).
"""
function _jacchia1977_temperature(z::Number, c::NTuple{7, T}) where {T <: Number}
    T₀ = _JACCHIA1977_CONSTANTS.T₀
    z₀ = _JACCHIA1977_CONSTANTS.z₀
    zx = _JACCHIA1977_CONSTANTS.zx

    Δzx = z - zx
    Δz₀ = z - z₀

    (Δz₀ == 0) && return T(T₀)

    if Δzx > 0
        return c[7] + c[4] * atan(c[5] * Δzx + c[6] * Δzx^3)
    else
        aux = Δzx / Δz₀
        return c[7] + c[1] * atan(c[2] * Δzx + c[3] * Δzx * aux * aux)
    end
end

"""
    _jacchia1977_static(T∞::T1, z::T2) where {T1<:Number, T2<:Number} -> NTuple{6, T}, T, T

Compute the Jacchia 1977 static model (routine IMOWEI of [2]) for the exospheric
temperature `T∞` [K] and altitude `z` [km].

The function numerically integrates the barometric equation between 90 km and 100 km and
the diffusion equations above 100 km using the Boole rule, as in the reference
implementation [2]. The atomic hydrogen is anchored at 500 km and integrated with its flux
term for other altitudes.

# Returns

- `NTuple{6, T}`: Base-10 logarithm of the number densities [1 / m³] in the internal order
    (He, O₂, N₂, Ar, O, H), where `T` is the promotion of `T1` and `T2`.
- `T`: Mean molecular mass at the selected altitude [g / mol].
- `T`: Total density at the selected altitude [kg / m³].

# References

- **[1]** Jacchia, L. G (1977). *Thermospheric temperature, density and composition: New
    models*. SAO Special Report #375.
- **[2]** de Matos, B. S., Carrara, V (1985-1987). *Fortran implementation of the Jacchia
    1977 model*. INPE, São José dos Campos, BR.
"""
function _jacchia1977_static(T∞::T1, z::T2) where {T1 <: Number, T2 <: Number}
    Rstar = _JACCHIA1977_CONSTANTS.Rstar
    Av    = _JACCHIA1977_CONSTANTS.Av
    Ra    = _JACCHIA1977_CONSTANTS.Ra
    g₀    = _JACCHIA1977_CONSTANTS.g₀
    M₀    = _JACCHIA1977_CONSTANTS.M₀
    qi    = _JACCHIA1977_CONSTANTS.qi
    T₀    = _JACCHIA1977_CONSTANTS.T₀
    ρ₀    = _JACCHIA1977_CONSTANTS.ρ₀
    Mi    = _JACCHIA1977_CONSTANTS.Mi
    αi    = _JACCHIA1977_CONSTANTS.αi
    Ca    = _JACCHIA1977_CONSTANTS.Ca
    Wb    = _JACCHIA1977_CONSTANTS.Wb

    RT = float(promote_type(T1, T2))

    ln10 = log(RT(10))

    c = _jacchia1977_profile_params(T∞)

    # Hydrogen flux and number density at 500 km (mks), Section 7 [1].
    aux       = 28.9 / T∞^RT(0.25)
    ϕH        = 10^(RT(6.90) + aux) / 2.0e20
    ln_nH_500 = (RT(5.94) + aux) * ln10

    # Number of species included in the mean molecular mass (H is included above 140 km).
    nc = 5

    # `an` contains the natural logarithm of the number densities in the internal order
    # (He, O₂, N₂, Ar, O, H).
    an = ntuple(_ -> RT(0), Val(6))

    ########################################################################################
    #           Barometric Equation Between 90 km and min(z, 100 km), Eq. 8 [1]            #
    ########################################################################################

    z_end = min(z, RT(100))

    int = zero(RT)
    zᵢ  = RT(90)

    step = _jacchia1977_step(zᵢ, z_end, RT(0.05))

    while abs(z_end - zᵢ) > 1e-4
        g  = g₀ / (1 + (zᵢ + 2step) / Ra)^2
        Σ  = zero(RT)
        zⱼ = zᵢ

        @inbounds for i in 1:5
            Δz = zⱼ - 90
            M′ = @evalpoly(Δz, Ca[1], Ca[2], Ca[3], Ca[4], Ca[5], Ca[6])
            Σ += Wb[i] * g * M′ / _jacchia1977_temperature(zⱼ, c)
            zⱼ += step
        end

        int += step * Σ
        zᵢ  += 4step
    end

    ρ′  = ρ₀ * exp(-int / Rstar)
    Δz  = z_end - 90
    M′  = @evalpoly(Δz, Ca[1], Ca[2], Ca[3], Ca[4], Ca[5], Ca[6])
    Tf  = _jacchia1977_temperature(z_end, c)
    N′  = Av * ρ′ / M₀ * T₀ / Tf
    ρ′  = N′ * M′
    aux = ρ′ / M₀

    @reset an[1] = log(qi[1] * aux)
    @reset an[2] = log(aux * (1 + qi[2]) - N′)
    @reset an[3] = log(qi[3] * aux)
    @reset an[4] = log(qi[4] * aux)
    @reset an[5] = log(2 * (N′ - aux))
    @reset an[6] = RT(0)

    if z > 100
        ####################################################################################
        #        Diffusion Equations Between 100 km and min(z, 140 km), Eq. 16 [1]         #
        ####################################################################################

        z_ini = RT(100)
        z_end = min(z, RT(140))
        Tᵢ    = Tf
        int   = zero(RT)

        step = _jacchia1977_step(z_ini, z_end, RT(0.05))
        zᵢ   = z_ini

        while abs(z_end - zᵢ) > 1e-4
            g  = g₀ / (1 + (zᵢ + 2step) / Ra)^2 / Rstar
            Σ  = zero(RT)
            zⱼ = zᵢ

            @inbounds for i in 1:5
                Σ  += Wb[i] * g / _jacchia1977_temperature(zⱼ, c)
                zⱼ += step
            end

            int += step * Σ
            zᵢ  += 4step
        end

        Tf  = _jacchia1977_temperature(z_end, c)
        aux = log(Tᵢ / Tf)

        for i in 1:5
            @reset an[i] = an[i] - int * Mi[i] + aux * (1 + αi[i])
        end

        if z > 140
            ################################################################################
            #     Diffusion Equations Between 140 km and 500 km (H Anchor), Eq. 16 [1]     #
            ################################################################################

            nc    = 6
            z_ini = RT(140)
            z_end = RT(500)
            Tᵢ    = Tf
            int   = zero(RT)
            step  = RT(5)
            zᵢ    = z_ini

            while abs(z_end - zᵢ) > 1e-4
                g  = g₀ / (1 + (zᵢ + 2step) / Ra)^2 / Rstar
                Σ  = zero(RT)
                zⱼ = zᵢ

                @inbounds for i in 1:5
                    Σ  += Wb[i] * g / _jacchia1977_temperature(zⱼ, c)
                    zⱼ += step
                end

                int += step * Σ
                zᵢ  += 4step
            end

            Tf  = _jacchia1977_temperature(z_end, c)
            aux = log(Tᵢ / Tf)

            @reset an[6] = ln_nH_500

            for i in 1:5
                @reset an[i] = an[i] - int * Mi[i] + aux * (1 + αi[i])
            end

            if z < 500
                ############################################################################
                #        Downward Integration From 500 km to z With the H Flux Term        #
                ############################################################################

                Tᵢ    = Tf
                z_ini = RT(500)
                z_end = RT(z)

                step = _jacchia1977_step(z_ini, z_end, RT(-5))
                zᵢ   = z_ini

                al = ntuple(_ -> RT(0), Val(5))

                while abs(z_end - zᵢ) > 1e-4
                    # Number densities with the departures from diffusive equilibrium
                    # corrections (eqs. 14 and 15 of [1]) at the beginning of the segment.
                    @reset al[1] = an[1]
                    @reset al[2] =
                        an[2] - RT(0.07) * (1 + tanh(RT(0.18) * (zᵢ - 111))) * ln10
                    @reset al[3] = an[3]
                    @reset al[4] = an[4]
                    @reset al[5] =
                        an[5] - RT(0.24) * exp(-RT(0.009) * (zᵢ - RT(97.7))^2) * ln10

                    g  = g₀ / (1 + (zᵢ + 2step) / Ra)^2 / Rstar
                    Σ  = zero(RT)
                    Σn = zero(RT)
                    zⱼ = zᵢ

                    @inbounds for i in 1:5
                        Σ  += Wb[i] * g / _jacchia1977_temperature(zⱼ, c)
                        Σn += exp(al[i])
                        zⱼ += step
                    end

                    zᵢ += 4step

                    int = step * Σ
                    Tf = _jacchia1977_temperature(zᵢ, c)
                    ΔlnT = log(Tᵢ / Tf)
                    Tᵢ = Tf

                    for i in 1:5
                        @reset an[i] = an[i] - int * Mi[i] + ΔlnT * (1 + αi[i])
                    end

                    # Hydrogen integration with the flux term, following the reference
                    # implementation [2].
                    Σϕ = Σn / exp(an[6]) * ϕH
                    g  = g₀ / (1 + (zᵢ + 2step) / Ra)^2 / Rstar
                    zⱼ = zᵢ
                    Σ  = zero(RT)
                    Σ₁ = zero(RT)
                    Σ₂ = zero(RT)

                    @inbounds for i in 1:5
                        Tl  = _jacchia1977_temperature(zⱼ, c)
                        Σ   += Wb[i] / Tl
                        aux = Wb[i] / √Tl
                        Σ₁  += Σϕ * aux
                        Σ₂  += aux
                        zⱼ  += step
                    end

                    @reset an[6] =
                        an[6] - (g * Σ * step * Mi[6] - ΔlnT * (1 + αi[6])) -
                        Σ₁ * 1000 * step - Σ₂ * ϕH * 1000 * step
                end
            elseif z > 500
                ############################################################################
                #         Upward Integration From 500 km to z With the H Flux Term         #
                ############################################################################

                Tᵢ    = Tf
                z_ini = RT(500)
                z_end = RT(z)
                int   = zero(RT)
                intϕ  = zero(RT)

                step = _jacchia1977_step(z_ini, z_end, RT(2.5))
                zᵢ   = z_ini

                while abs(z_end - zᵢ) > 1e-4
                    g  = g₀ / (1 + (zᵢ + 2step) / Ra)^2 / Rstar
                    Σ  = zero(RT)
                    Σ₂ = zero(RT)
                    zⱼ = zᵢ

                    @inbounds for i in 1:5
                        Tl = _jacchia1977_temperature(zⱼ, c)
                        Σ += Wb[i] * g / Tl
                        Σ₂ += Wb[i] / √Tl
                        zⱼ += step
                    end

                    int  += step * Σ
                    intϕ += 1000 * step * ϕH * Σ₂
                    zᵢ   += 4step
                end

                Tf  = _jacchia1977_temperature(z_end, c)
                aux = log(Tᵢ / Tf)

                for i in 1:6
                    @reset an[i] = an[i] - int * Mi[i] + aux * (1 + αi[i])
                end

                @reset an[6] = an[6] - intϕ
            end
        end
    end

    ########################################################################################
    #        Departures From Diffusive Equilibrium (Eqs. 14 and 15 [1]) and Output         #
    ########################################################################################

    @reset an[2] = an[2] - RT(0.07) * (1 + tanh(RT(0.18) * (z - 111))) * ln10
    @reset an[5] = an[5] - RT(0.24) * exp(-RT(0.009) * (z - RT(97.7))^2) * ln10

    Σm = zero(RT)
    Σn = zero(RT)
    log₁₀_n = ntuple(_ -> RT(0), Val(6))

    @inbounds for i in 1:6
        if i <= nc
            nᵢ = exp(an[i])
            Σm += nᵢ * Mi[i]
            Σn += nᵢ
        end

        # Convert to base-10 logarithm and clamp negative values, as in the reference
        # implementation [2].
        @reset log₁₀_n[i] = max(an[i] / ln10, RT(0))
    end

    M̄ = Σm / Σn
    ρ = Σm / Av

    return log₁₀_n, M̄, ρ
end

"""
    _jacchia1977_step(z_ini::Number, z_end::Number, base_step::Number) -> Number

Compute the integration step [km] so that the interval between `z_ini` and `z_end` [km] is
divided into an integer number of Boole rule applications with a step close to `base_step`
[km], as in the reference implementation [2].
"""
function _jacchia1977_step(z_ini::Number, z_end::Number, base_step::Number)
    quarter = (z_end - z_ini) / 4
    n = trunc(quarter / base_step)
    (n <= 0) && (n = one(n))
    return quarter / n
end

"""
    _jacchia1977_diurnal(
        T½::Number,
        Ωp::Number,
        Ωs::Number,
        δs::Number,
        ϕ::Number,
        z::Number,
        M̄::Number
    ) -> NTuple{6, T}, NTuple{6, T}, T

Compute the base-10 logarithm of the number densities considering the diurnal variation
(routine DIVARI of [2] and eqs. 24 to 27 of [1]).

# Arguments

- `T½::Number`: Mean exospheric temperature [K].
- `Ωp::Number`: Right ascension of the selected location [rad].
- `Ωs::Number`: Right ascension of the Sun [rad].
- `δs::Number`: Declination of the Sun [rad].
- `ϕ::Number`: Latitude [rad].
- `z::Number`: Altitude [km].
- `M̄::Number`: Local mean molecular mass [g / mol].

# Returns

- `NTuple{6, T}`: Base-10 logarithm of the number densities [1 / m³] evaluated at the
    per-species pseudo exospheric temperatures, in the internal order (He, O₂, N₂, Ar, O,
    H), where `T` is the promotion of the input types.
- `NTuple{6, T}`: Base-10 logarithm of the number densities [1 / m³] of the static model
    evaluated at the hydrogen pseudo exospheric temperature.
- `T`: Pseudo exospheric temperature of the hydrogen [K], used as the quiet temperature by
    the geomagnetic variation.
"""
function _jacchia1977_diurnal(
    T½::Number, Ωp::Number, Ωs::Number, δs::Number, ϕ::Number, z::Number, M̄::Number
)
    Mi = _JACCHIA1977_CONSTANTS.Mi

    RT = float(
        promote_type(
            typeof(T½), typeof(Ωp), typeof(Ωs), typeof(δs), typeof(ϕ), typeof(z), typeof(M̄)
        ),
    )

    # Hour angle of the Sun [rad].
    H = Ωp - Ωs

    # Exponent of the diurnal variation, eq. 26 [1].
    n = 2 + cos(ϕ * ϕ / (π / 2))^2

    aux = 1 + RT(0.3666069) * δs * sin(ϕ)
    cos_ϕ = cos(ϕ)

    al  = ntuple(_ -> RT(0), Val(6))
    ac  = ntuple(_ -> RT(0), Val(6))
    Θ_H = RT(0)

    @inbounds for i in 1:6
        # Phase angle of the diurnal variation, eq. 27 [1]. For the hydrogen, the phase
        # angle is -60°.
        β = (i == 6) ? -RT(1.0471976) : deg2rad(27 * (M̄ / Mi[i] - 1) - 35)

        A = H + β
        f = RT(0.08) * cos(3A - RT(1.3089969)) + abs(cos(A / 2))^n

        # Pseudo exospheric temperature of the i-th species, eq. 24 [1].
        Θᵢ = T½ * (aux + RT(0.24) * cos_ϕ * (f - RT(0.5)))

        an, _, _ = _jacchia1977_static(Θᵢ, z)

        @reset al[i] = an[i]

        if i == 6
            ac  = an
            Θ_H = Θᵢ
        end
    end

    return al, ac, Θ_H
end

"""
    _jacchia1977_geomagnetic(
        T_quiet::Number,
        Kp::Number,
        ϕ::Number,
        λ::Number,
        z::Number
    ) -> NTuple{6, T}

Compute the base-10 logarithm of the number densities of the static model evaluated at the
exospheric temperature increased by the geomagnetic activity `Kp` [-], together with the
homopause displacement and equatorial wave corrections (routine GEOACI of [2] and eqs. 28
to 35 of [1]), given the quiet exospheric temperature `T_quiet` [K], the latitude `ϕ`
[rad], the longitude `λ` [rad], and the altitude `z` [km].

# Returns

- `NTuple{6, T}`: Base-10 logarithm of the number densities [1 / m³] in the internal order
    (He, O₂, N₂, Ar, O, H), where `T` is the promotion of the input types. The caller must
    subtract the static model evaluated at `T_quiet` to obtain the geomagnetic variation.
"""
function _jacchia1977_geomagnetic(
    T_quiet::Number, Kp::Number, ϕ::Number, λ::Number, z::Number
)
    ai = _JACCHIA1977_CONSTANTS.ai

    # Amplitude of the geomagnetic effect, eq. 31 [1].
    A = 57.5 * Kp * (1 + 0.027 * exp(0.4 * Kp))

    # Sine of the geomagnetic (invariant) latitude using the dipole approximation [2].
    sin_ϕᵢ  = 0.9792 * sin(ϕ) + 0.2028 * cos(ϕ) * cos(λ - 5.0789081)
    sin²_ϕᵢ = sin_ϕᵢ * sin_ϕᵢ
    cos²_ϕᵢ = 1 - sin²_ϕᵢ

    # Geomagnetic variation of the exospheric temperature.
    ΔT∞ = A * sin²_ϕᵢ * sin²_ϕᵢ

    dn, _, _ = _jacchia1977_static(T_quiet + ΔT∞, z)

    # Homopause displacement, eq. 33 [1] [m].
    aux  = 0.01 * ΔT∞
    Δz_H = 5000 * log(√(1 + aux * aux) + aux)

    # Equatorial wave, eq. 35 [1].
    Δe = 5.2e-4 * A * cos²_ϕᵢ * cos²_ϕᵢ

    return ntuple(i -> dn[i] + ai[i] * Δz_H + Δe, Val(6))
end

"""
    _jacchia1977_seasonal_latitudinal(
        Φ::Number,
        δs::Number,
        ϕ::Number,
        z::Number
    ) -> NTuple{6, T}

Compute the seasonal-latitudinal variation of the base-10 logarithm of the number
densities (routine SEALAT of [2] and eqs. 36 to 39 of [1]) given the fraction of the
tropic year `Φ` [-], the Sun declination `δs` [rad], the latitude `ϕ` [rad], and the
altitude `z` [km]. The result is in the internal order (He, O₂, N₂, Ar, O, H), and `T` is
the promotion of the input types.
"""
function _jacchia1977_seasonal_latitudinal(Φ::Number, δs::Number, ϕ::Number, z::Number)
    ci = _JACCHIA1977_CONSTANTS.ci

    sin_ϕ = sin(ϕ)

    # Thermospheric variation, eq. 36 [1].
    Δt = δs * sin_ϕ / 0.409157536545

    # Mesospheric variation, eqs. 37 to 39 [1].
    #
    # Notice that the reference implementation [2] applies the sign of the latitude to the
    # absolute value of the entire term due to the Fortran DSIGN intrinsic, dropping the
    # sign of `P`. We follow eq. 37 of the report [1] instead, in which only the sign of
    # the latitude multiplies the term.
    Δz = z - 91
    S  = 0.014 * Δz * exp(-0.0013 * Δz * Δz)
    P  = sin(2π * Φ + 1.72)
    Δm = flipsign(sin_ϕ * sin_ϕ * S * P, ϕ)

    return ntuple(i -> Δt * ci[i] + Δm, Val(6))
end

"""
    _jacchia1977_semiannual(Φ::Number, z::Number) -> Number

Compute the semiannual variation of the base-10 logarithm of the number densities (routine
SEMIAN of [2] and eqs. 40 to 44 of [1]) given the fraction of the tropic year `Φ` [-] and
the altitude `z` [km]. The variation is the same for all species.
"""
function _jacchia1977_semiannual(Φ::Number, z::Number)
    f = (0.04 * z * z / 1e4 + 0.05) * exp(-0.0025 * z)
    τ = Φ + 0.0954 * ((0.5 + 0.5 * sin(2π * Φ + 6.04))^1.65 - 0.5)
    g = 0.0284 + 0.382 * sin(4π * τ + 4.26) * (1 + 0.467 * sin(2π * τ + 4.14))

    return f * g
end
