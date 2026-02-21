export harrispriester

"""
    harrispriester(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number; kwargs...) -> Float64
    harrispriester(jd::Number, ϕ_gd::Number, λ::Number, h::Number; kwargs...) -> Float64

Compute the atmospheric density [kg / m³] using the Harris-Priester model.

# Arguments
- `jd`: The Julian day to compute the model.
- `ϕ_gd`: Geodetic latitude [rad].
- `λ`: Geodetic longitude [rad].
- `h`: Geodetic altitude [m].
- `n`: Cosine exponent in the diurnal bulge modeling (`2` <= `n` <= `6`).
    - If `n` is `2`, it models a smooth transition.
    - If `n` is `6`, it models a sharp transition.
    - **Default**: 4
- `alt_ρ::AbstractMatrix`: A matrix containing the minimum and maximum density profiles.
    - `[:, 1]`: Altitude [m].
    - `[:, 2]`: Minimum density [kg / m³].
    - `[:, 3]`: Maximum density [kg / m³].
    - **Default**: `_HARRIS_PRIESTER_ALT_RHO` (mean solar activity).

# Returns

- The atmospheric density [kg / m³].
"""
function harrispriester(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number;
    n::Int = 4,
    alt_ρ::AbstractMatrix{<:Number} = _HARRIS_PRIESTER_ALT_RHO
)
    return harrispriester(datetime2julian(instant), ϕ_gd, λ, h; n = n, alt_ρ = alt_ρ)
end

function harrispriester(
    jd::JT,
    ϕ_gd::PT,
    λ::LT,
    h::HT;
    n::Int = 4,
    alt_ρ::AbstractMatrix{DT} = _HARRIS_PRIESTER_ALT_RHO
) where {JT<:Number, PT<:Number, LT<:Number, HT<:Number, DT<:Number}

    @assert 2 <= n <= 6 "The cosine exponent must be between 2 and 6."

    RT = promote_type(JT, PT, LT, HT, DT)

    min_alt = alt_ρ[1, 1]
    max_alt = alt_ρ[end, 1]

    if h < min_alt
        throw(ArgumentError("The altitude is lower than the minimum altitude in the density profile."))
    elseif h > max_alt
        return RT(0)
    end

    # Compute the Sun position represented in the inertial reference frame (MOD).
    s_i = sun_position_mod(jd)

    # Compute the Sun declination [rad].
    δs = atan(s_i[3], √(s_i[1]^2 + s_i[2]^2))

    # Compute the Sun right ascension [rad].
    Ωs = atan(s_i[2], s_i[1])

    # Compute the right ascension of the selected location w.r.t. the inertial reference
    # frame.
    Ωp = λ + jd_to_gmst(jd)

    # Compute the cosine of the angle between the diurnal bulge apex and the satellite.
    #
    # The bulge is lagged by 30 degrees (approx 2 hours) from the Sun's RA. We assume the
    # bulge has the same declination as the Sun.
    #
    #   cos(ψ) = sin(δ_sat)sin(δ_bulge) + cos(δ_sat)cos(δ_bulge)cos(RA_sat - RA_bulge)
    #
    #   δ_sat    = ϕ_gd
    #   RA_sat   = Ωp
    #   δ_bulge  = δs
    #   RA_bulge = Ωs + lag

    sin_δs, cos_δs = sincos(δs)
    sin_ϕ, cos_ϕ   = sincos(ϕ_gd)

    cos_ψ = sin_ϕ * sin_δs + cos_ϕ * cos_δs * cos(Ωp - Ωs - _HARRIS_PRIESTER_LAG_ANGLE)
    c2ψ2 = (1 + cos_ψ) / 2
    cψ2 = √(c2ψ2)
    cos_pow = cψ2 > _HARRIS_PRIESTER_MIN_COS ? cψ2^n : RT(0)

    # Search for the altitude index in the density table.
    alt_col = @view alt_ρ[:, 1]
    ia = searchsortedlast(alt_col, h)
    ia = max(1, min(ia, size(alt_ρ, 1) - 1))

    # Fractional satellite height.
    h₁ = alt_ρ[ia, 1]
    h₂ = alt_ρ[ia + 1, 1]
    dh = (h₁ - h) / (h₁ - h₂)

    # Minimum exponential density interpolation.
    ρ_min₁ = alt_ρ[ia, 2]
    ρ_min₂ = alt_ρ[ia + 1, 2]
    ρ_min = ρ_min₁ * (ρ_min₂ / ρ_min₁)^dh

    if abs(cos_pow) < eps(RT)
        return ρ_min
    else
        # Maximum exponential density interpolation.
        ρ_max₁ = alt_ρ[ia, 3]
        ρ_max₂ = alt_ρ[ia + 1, 3]
        ρ_max = ρ_max₁ * (ρ_max₂ / ρ_max₁)^dh
        return ρ_min + (ρ_max - ρ_min) * cos_pow
    end
end
