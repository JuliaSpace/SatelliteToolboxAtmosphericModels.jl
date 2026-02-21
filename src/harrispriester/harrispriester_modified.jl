export harrispriester_modified

"""
    harrispriester_modified(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number, F10ₐ::Number; n::Number = 4) -> Float64
    harrispriester_modified(jd::Number, ϕ_gd::Number, λ::Number, h::Number, F10ₐ::Number; n::Number = 4) -> Float64
    harrispriester_modified(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number; n::Number = 4) -> Float64
    harrispriester_modified(jd::Number, ϕ_gd::Number, λ::Number, h::Number; n::Number = 4) -> Float64

Compute the atmospheric density [kg / m³] using the modified Harris-Priester model.

This model is a Julia translation of the Fortran code `harris_priester_mod_dist.f90`
developed by Noble Hatten and Ryan P. Russell [1]. It ensures continuous first derivatives,
eliminates singularities, and uses a cubic dependency on the 81-day centered average of the
F10.7 solar flux index (`F10ₐ`) to model density variations.

If `F10ₐ` is not provided, it will be automatically fetched using the `SpaceIndices` package.
In this case, the initialization of the space indices package with `SpaceIndices.init()` is
required.

# References

- **[1]** Hatten, N., & Russell, R. P. (2017). A smooth and robust Harris-Priester
    atmospheric density model for low Earth orbit applications. *Advances in Space
    Research*, 59(2), 571-586.

# Arguments

- `instant`: The `DateTime` to compute the model.
- `jd`: The Julian day to compute the model.
- `ϕ_gd`: Geodetic latitude [rad].
- `λ`: Geodetic longitude [rad].
- `h`: Geodetic altitude [m].
- `F10ₐ`: (Optional) 81-day centered average of the F10.7 solar flux index [sfu].
- `n`: (Optional) Cosine exponent in the diurnal bulge modeling. Default: 4.
    **IMPORTANT:** The original FORTRAN implementation computes `n` from orbital inclination:
    ```
    n = 2.001 + 4 * sin²(inclination)
    ```
    - For equatorial orbits (inclination = 0°): use `n = 2.001`
    - For polar orbits (inclination = 90°): use `n = 6.001`
    - For intermediate orbits: calculate appropriately
    If orbital inclination is unknown, the default `n = 4` provides a reasonable approximation. This functionality was
    purposefully removed here to avoid needing an additional dependency.

# Returns

- The atmospheric density [kg / m³].
"""
function harrispriester_modified(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number,
    F10ₐ::Number;
    n::Number = 4
)
    return harrispriester_modified(datetime2julian(instant), ϕ_gd, λ, h, F10ₐ; n = n)
end

function harrispriester_modified(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number;
    n::Number = 4
)
    return harrispriester_modified(datetime2julian(instant), ϕ_gd, λ, h; n = n)
end

function harrispriester_modified(
    jd::Number,
    ϕ_gd::Number,
    λ::Number,
    h::Number;
    n::Number = 4
)
    # Fetch the 81-day centered average of F10.7 solar flux.
    F10ₐ = sum(space_index.(Val(:F10obs), jd + k) for k in -40:40) / 81
    return harrispriester_modified(jd, ϕ_gd, λ, h, F10ₐ; n = n)
end

function harrispriester_modified(
    jd::JT,
    ϕ_gd::PT,
    λ::LT,
    h::HT,
    F10ₐ::FT;
    n::Number = 4
) where {JT<:Number, PT<:Number, LT<:Number, HT<:Number, FT<:Number}

    RT = promote_type(JT, PT, LT, HT, FT, typeof(n))

    # Convert inputs to kilometers for consistency with the original Fortran model.
    h_km = h / 1000

    # Get altitude index.
    if h_km < _HARRIS_PRIESTER_MOD_HVEC[1]
        i = 1
    elseif h_km >= _HARRIS_PRIESTER_MOD_HVEC[end]
        i = length(_HARRIS_PRIESTER_MOD_HVEC)
    else
        i = searchsortedlast(_HARRIS_PRIESTER_MOD_HVEC, h_km)
    end
    
    if i < length(_HARRIS_PRIESTER_MOD_HVEC)
        ρ_maxᵢ   = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i, 1:4])...)
        ρ_minᵢ   = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i, 5:8])...)
        ρ_maxᵢ₊₁ = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i + 1, 1:4])...)
        ρ_minᵢ₊₁ = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i + 1, 5:8])...)
        hᵢ         = _HARRIS_PRIESTER_MOD_HVEC[i]
        hᵢ₊₁       = _HARRIS_PRIESTER_MOD_HVEC[i + 1]
    else # Above nominal altitude max
        ρ_maxᵢ   = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i - 1, 1:4])...)
        ρ_minᵢ   = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i - 1, 5:8])...)
        ρ_maxᵢ₊₁ = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i, 1:4])...)
        ρ_minᵢ₊₁ = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i, 5:8])...)
        hᵢ         = _HARRIS_PRIESTER_MOD_HVEC[i-1]
        hᵢ₊₁       = _HARRIS_PRIESTER_MOD_HVEC[i]
    end

    Δhᵢ = hᵢ₊₁ - hᵢ

    # Scale heights are calculated by exponential interpolation to maintain continuity in
    # density between altitude layers [1, Eq. 4, 5].
    H_ρ_minᵢ = -Δhᵢ / log(ρ_minᵢ₊₁ / ρ_minᵢ)
    H_ρ_maxᵢ = -Δhᵢ / log(ρ_maxᵢ₊₁ / ρ_maxᵢ)

    # A polynomial weighting function is used to ensure third-order continuity of scale
    # heights when passing through altitude boundaries [1, Section 3.2].
    α = _HARRIS_PRIESTER_MOD_α

    if (h_km <= hᵢ + α) && (i > 1) # Near lower boundary
        hᵢ₋₁       = _HARRIS_PRIESTER_MOD_HVEC[i-1]
        ρ_maxᵢ₋₁ = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i - 1, 1:4])...)
        ρ_minᵢ₋₁ = @evalpoly(F10ₐ, @view(_HARRIS_PRIESTER_MOD_COEFS[i - 1, 5:8])...)

        Δhᵢ₋₁ = hᵢ - hᵢ₋₁
        xbar     = SVector{2}(hᵢ - α, hᵢ + α)

        H_ρ_minᵢ₋₁ = -Δhᵢ₋₁ / log(ρ_minᵢ / ρ_minᵢ₋₁)
        H_ρ_maxᵢ₋₁ = -Δhᵢ₋₁ / log(ρ_maxᵢ / ρ_maxᵢ₋₁)

        H_ρ_min_vec = SVector{2}(H_ρ_minᵢ₋₁, H_ρ_minᵢ)
        H_ρ_max_vec = SVector{2}(H_ρ_maxᵢ₋₁, H_ρ_maxᵢ)

        _, H_ρ_min′ᵢ, H_ρ_max′ᵢ =
            _scale_height_junk(xbar, h_km, H_ρ_min_vec, H_ρ_max_vec)

        ρ_min_h = ρ_minᵢ * exp((hᵢ - h_km) / H_ρ_min′ᵢ)
        ρ_max_h = ρ_maxᵢ * exp((hᵢ - h_km) / H_ρ_max′ᵢ)

    elseif (h_km >= hᵢ₊₁ - α) && (i < length(_HARRIS_PRIESTER_MOD_HVEC) - 1) # Near upper boundary
        hᵢ₊₂       = _HARRIS_PRIESTER_MOD_HVEC[i+2]
        ρ_maxᵢ₊₂ = @evalpoly(F10ₐ, view(_HARRIS_PRIESTER_MOD_COEFS, i + 2, 1:4)...)
        ρ_minᵢ₊₂ = @evalpoly(F10ₐ, view(_HARRIS_PRIESTER_MOD_COEFS, i + 2, 5:8)...)

        Δhᵢ₊₁ = hᵢ₊₂ - hᵢ₊₁
        xbar       = SVector{2}(hᵢ₊₁ - α, hᵢ₊₁ + α)

        H_ρ_minᵢ₊₁ = -Δhᵢ₊₁ / log(ρ_minᵢ₊₂ / ρ_minᵢ₊₁)
        H_ρ_maxᵢ₊₁ = -Δhᵢ₊₁ / log(ρ_maxᵢ₊₂ / ρ_maxᵢ₊₁)

        H_ρ_min_vec = SVector{2}(H_ρ_minᵢ, H_ρ_minᵢ₊₁)
        H_ρ_max_vec = SVector{2}(H_ρ_maxᵢ, H_ρ_maxᵢ₊₁)

        _, H_ρ_min′ᵢ, H_ρ_max′ᵢ =
            _scale_height_junk(xbar, h_km, H_ρ_min_vec, H_ρ_max_vec)

        ρ_min_h = ρ_minᵢ₊₁ * exp((hᵢ₊₁ - h_km) / H_ρ_min′ᵢ)
        ρ_max_h = ρ_maxᵢ₊₁ * exp((hᵢ₊₁ - h_km) / H_ρ_max′ᵢ)

    else # Not near a boundary, no weighting needed.
        ρ_min_h = ρ_minᵢ * exp((hᵢ - h_km) / H_ρ_minᵢ)
        ρ_max_h = ρ_maxᵢ * exp((hᵢ - h_km) / H_ρ_maxᵢ)
    end

    # Compute the Sun position represented in the inertial reference frame (MOD).
    s_i = sun_position_mod(jd)

    # Compute the Sun declination [rad].
    δs = atan(s_i[3], √(s_i[1]^2 + s_i[2]^2))

    # Compute the Sun right ascension [rad].
    Ωs = atan(s_i[2], s_i[1])

    # Compute the right ascension of the selected location w.r.t. the inertial reference frame.
    Ωp = λ + jd_to_gmst(jd)

    # Compute the cosine of the angle between the diurnal bulge apex and the satellite.
    #
    # The bulge is lagged by 30 degrees (approx 2 hours) from the Sun's RA. We assume the
    # bulge has the same declination as the Sun.
    #
    #   cos(ψ) = sin(δ_sat)sin(δ_bulge) + cos(δ_sat)cos(δ_bulge)cos(RA_sat - RA_bulge)

    sin_δs, cos_δs = sincos(δs)
    sin_ϕ, cos_ϕ   = sincos(ϕ_gd)

    cos_ψ = sin_ϕ * sin_δs + cos_ϕ * cos_δs * cos(Ωp - Ωs - _HARRIS_PRIESTER_LAG_ANGLE)
    cos²_ψ_by_2 = max(RT(0), 0.5 + 0.5 * cos_ψ)
    cos_ψ_by_2 = sqrt(cos²_ψ_by_2)

    # Compute cos(ψ/2)ⁿ with a smoothing function to avoid issues when cos(ψ/2) is near 0.
    cos_ψ_by_2_pow_n = RT(0)
    if cos_ψ_by_2 >= _HARRIS_PRIESTER_MOD_COS_ψ_BY_2_TOL
        cos_ψ_by_2_pow_n = cos_ψ_by_2^n
    elseif cos_ψ_by_2 > 0
        ξ = cos_ψ_by_2 / _HARRIS_PRIESTER_MOD_COS_ψ_BY_2_TOL
        xiterm = @evalpoly(ξ, 10, -15, 6)
        c1 = cos_ψ_by_2^n
        w1 = ξ^3 * xiterm
        cos_ψ_by_2_pow_n = w1 * c1
    end

    # Final density calculation [1, Eq. 1].
    ρ = ρ_min_h + (ρ_max_h - ρ_min_h) * cos_ψ_by_2_pow_n

    # Convert from g/km³ to kg/m³.
    return ρ * 1e-12
end

# Private function to calculate scale heights via Junkins/Jancaitis weighting method for
# third-order continuity [1, Section 3.2].
function _scale_height_junk(
    xbar::AbstractVector{<:Number},
    h::Number,
    H_ρ_min::AbstractVector{<:Number},
    H_ρ_max::AbstractVector{<:Number}
)
    xbardiff = xbar[2] - xbar[1]
    ξ = (h - xbar[1]) / xbardiff

    # Weighting function for third-order continuity [1, Eq. 20].
    t1 = @evalpoly(ξ, 35, -84, 70, -20)
    w1 = ξ^4 * t1

    H_ρ_min′ = H_ρ_min[1] + w1 * (H_ρ_min[2] - H_ρ_min[1])
    H_ρ_max′ = H_ρ_max[1] + w1 * (H_ρ_max[2] - H_ρ_max[1])

    return w1, H_ρ_min′, H_ρ_max′
end
