## Description #############################################################################
#
# Modified Harris-Priester atmospheric density model.
#
############################################################################################

export harrispriester_modified

"""
    harrispriester_modified(
        instant::DateTime,
        ϕ_gd::Number,
        λ::Number,
        h::Number[, F10ₐ::Number];
        kwargs...
    ) -> Number
    harrispriester_modified(
        jd::Number,
        ϕ_gd::Number,
        λ::Number,
        h::Number[, F10ₐ::Number];
        kwargs...
    ) -> Number

Compute the atmospheric density [kg / m³] using the modified Harris-Priester model.

This model is a Julia translation of the Fortran code `harris_priester_mod_dist.f90`
developed by Noble Hatten and Ryan P. Russell [1]. It ensures continuous first derivatives,
eliminates singularities, and uses a cubic dependency on the 81-day centered average of the
F10.7 solar flux index (`F10ₐ`) to model density variations.

If `F10ₐ` is not provided, it will be automatically fetched using the `SpaceIndices`
package. In this case, the initialization of the space indices package with
`SpaceIndices.init()` is required.

The model is valid only between 100 km and 1000 km. The function throws an `ArgumentError`
if the altitude `h` is outside this range.

# Arguments

- `instant::DateTime`: Instant to compute the model represented using `DateTime`.
- `jd::Number`: Julian day to compute the model.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Geodetic longitude [rad].
- `h::Number`: Geodetic altitude [m].
- `F10ₐ::Number`: 81-day centered average of the F10.7 solar flux index [sfu].

# Keywords

- `n::Number`: Cosine exponent in the diurnal bulge modeling (`2 <= n <= 7`). The original
    Fortran implementation computes `n` from the orbital inclination as
    `n = 2.001 + 4 sin²(inclination)`. Hence, use `n = 2.001` for equatorial orbits,
    `n = 6.001` for polar orbits, and interpolate for intermediate inclinations. If the
    orbital inclination is unknown, the default provides a reasonable approximation. This
    functionality was purposefully removed here to avoid needing an additional dependency.
    The function throws an `ArgumentError` if `n` is outside this range.
    (**Default**: `4`)

# Returns

- `Number`: Atmospheric density [kg / m³]. The type is the promotion of the types of the
    numeric inputs.

# References

- **[1]** Hatten, N., & Russell, R. P. (2017). A smooth and robust Harris-Priester
    atmospheric density model for low Earth orbit applications. *Advances in Space
    Research*, 59(2), 571-586.
"""
function harrispriester_modified(
    instant::DateTime, ϕ_gd::Number, λ::Number, h::Number, F10ₐ::Number; n::Number = 4
)
    return harrispriester_modified(datetime2julian(instant), ϕ_gd, λ, h, F10ₐ; n = n)
end

function harrispriester_modified(
    instant::DateTime, ϕ_gd::Number, λ::Number, h::Number; n::Number = 4
)
    return harrispriester_modified(datetime2julian(instant), ϕ_gd, λ, h; n = n)
end

function harrispriester_modified(
    jd::Number, ϕ_gd::Number, λ::Number, h::Number; n::Number = 4
)
    _check_altitude(h, _HARRIS_PRIESTER_MOD_H_MIN, _HARRIS_PRIESTER_MOD_H_MAX)

    # Fetch the 81-day centered average of F10.7 solar flux.
    #
    # NOTE: The reference [1] does not state whether the fit uses the observed flux or the
    # flux adjusted to 1 AU. We use the observed flux, which is the series usually
    # provided to the Harris-Priester model in orbit propagators.
    F10ₐ = _f10_81day_mean(Val(:F10obs), jd)

    @debug """
    Modified Harris-Priester - Fetched Space Indices
      81-day averaged F10.7 : $(F10ₐ) sfu
    """

    return harrispriester_modified(jd, ϕ_gd, λ, h, F10ₐ; n = n)
end

function harrispriester_modified(
    jd::JT, ϕ_gd::PT, λ::LT, h::HT, F10ₐ::FT; n::Number = 4
) where {JT <: Number, PT <: Number, LT <: Number, HT <: Number, FT <: Number}
    (2 <= n <= 7) || throw(ArgumentError("The cosine exponent must be between 2 and 7."))

    _check_altitude(h, _HARRIS_PRIESTER_MOD_H_MIN, _HARRIS_PRIESTER_MOD_H_MAX)

    RT = float(promote_type(JT, PT, LT, HT, FT))

    # Convert inputs to kilometers for consistency with the original Fortran model.
    h_km = h / 1000

    # Get the altitude layer index. At the top of the table, we use the last layer.
    i = searchsortedlast(_HARRIS_PRIESTER_MOD_HVEC, h_km)
    i = min(i, length(_HARRIS_PRIESTER_MOD_HVEC) - 1)

    # == Diurnal Variation =================================================================

    # Compute the Sun declination, the Sun right ascension, and the right ascension of the
    # selected location [rad].
    δs, Ωs, Ωp = _sun_geometry(jd, λ)

    # Compute the cosine of the angle between the diurnal bulge apex and the satellite.
    #
    # The bulge is lagged by 30 degrees (approx 2 hours) from the Sun's RA. We assume the
    # bulge has the same declination as the Sun.
    #
    #   cos(ψ) = sin(δ_sat)sin(δ_bulge) + cos(δ_sat)cos(δ_bulge)cos(RA_sat - RA_bulge)

    sin_δs, cos_δs = sincos(δs)
    sin_ϕ, cos_ϕ   = sincos(ϕ_gd)

    cos_ψ = sin_ϕ * sin_δs + cos_ϕ * cos_δs * cos(Ωp - Ωs - _HARRIS_PRIESTER_LAG_ANGLE)
    cos²_ψ_by_2 = max(RT(0), RT(1 / 2) + RT(1 / 2) * cos_ψ)
    cos_ψ_by_2 = sqrt(cos²_ψ_by_2)

    # Compute cos(ψ/2)ⁿ with a smoothing function to avoid issues when cos(ψ/2) is near 0.
    cos_ψ_by_2_pow_n = RT(0)
    if cos_ψ_by_2 >= _HARRIS_PRIESTER_MOD_COS_ψ_BY_2_TOL
        cos_ψ_by_2_pow_n = RT(cos_ψ_by_2^n)
    elseif cos_ψ_by_2 > 0
        ξ = cos_ψ_by_2 / _HARRIS_PRIESTER_MOD_COS_ψ_BY_2_TOL
        xiterm = @evalpoly(ξ, 10, -15, 6)
        c1 = cos_ψ_by_2^n
        w1 = ξ^3 * xiterm
        cos_ψ_by_2_pow_n = RT(w1 * c1)
    end

    # == Density ===========================================================================

    # Minimum density profile at the altitude. On the night side, the diurnal factor is 0
    # and the maximum density profile is not required.
    ρ_min_h = _harris_priester_mod_profile_density(F10ₐ, h_km, i, 5)

    iszero(cos_ψ_by_2_pow_n) && return RT(ρ_min_h * 1e-12)

    ρ_max_h = _harris_priester_mod_profile_density(F10ₐ, h_km, i, 1)

    # Final density calculation [1, Eq. 1].
    ρ = ρ_min_h + (ρ_max_h - ρ_min_h) * cos_ψ_by_2_pow_n

    # Convert from g/km³ to kg/m³.
    return RT(ρ * 1e-12)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _harris_priester_mod_profile_density(
        F10ₐ::Number,
        h_km::Number,
        i::Integer,
        c::Integer
    ) -> Number

Compute the density [g / km³] of the minimum (`c = 5`) or maximum (`c = 1`) density profile
at the altitude `h_km` [km] inside the layer `i` of `_HARRIS_PRIESTER_MOD_HVEC` for the
81-day averaged F10.7 flux `F10ₐ` [sfu], where `c` is the first column of the polynomial
coefficients of the profile in `_HARRIS_PRIESTER_MOD_COEFS`.

The scale heights are obtained by exponential interpolation between the layer boundaries
[1, eqs. 4 and 5], and they are blended with the scale heights of the adjacent layers near
the boundaries using the Junkins / Jancaitis weighting for third-order continuity [1,
Section 3.2].
"""
function _harris_priester_mod_profile_density(
    F10ₐ::Number, h_km::Number, i::Integer, c::Integer
)
    HVEC = _HARRIS_PRIESTER_MOD_HVEC
    α    = _HARRIS_PRIESTER_MOD_α

    ρᵢ   = _harris_priester_mod_density_poly(F10ₐ, i, c)
    ρᵢ₊₁ = _harris_priester_mod_density_poly(F10ₐ, i + 1, c)
    hᵢ   = HVEC[i]
    hᵢ₊₁ = HVEC[i + 1]

    # Scale height of the layer by exponential interpolation [1, Eq. 4, 5].
    Hᵢ = -(hᵢ₊₁ - hᵢ) / log(ρᵢ₊₁ / ρᵢ)

    if (h_km <= hᵢ + α) && (i > 1)
        # Near the lower boundary, blend with the scale height of the previous layer.
        ρᵢ₋₁ = _harris_priester_mod_density_poly(F10ₐ, i - 1, c)
        hᵢ₋₁ = HVEC[i - 1]
        Hᵢ₋₁ = -(hᵢ - hᵢ₋₁) / log(ρᵢ / ρᵢ₋₁)
        H′   = _scale_height_junk(hᵢ - α, hᵢ + α, h_km, Hᵢ₋₁, Hᵢ)

        return ρᵢ * exp((hᵢ - h_km) / H′)

    elseif (h_km >= hᵢ₊₁ - α) && (i < length(HVEC) - 1)
        # Near the upper boundary, blend with the scale height of the next layer.
        ρᵢ₊₂ = _harris_priester_mod_density_poly(F10ₐ, i + 2, c)
        hᵢ₊₂ = HVEC[i + 2]
        Hᵢ₊₁ = -(hᵢ₊₂ - hᵢ₊₁) / log(ρᵢ₊₂ / ρᵢ₊₁)
        H′   = _scale_height_junk(hᵢ₊₁ - α, hᵢ₊₁ + α, h_km, Hᵢ, Hᵢ₊₁)

        return ρᵢ₊₁ * exp((hᵢ₊₁ - h_km) / H′)
    end

    # Not near a boundary, no blending needed.
    return ρᵢ * exp((hᵢ - h_km) / Hᵢ)
end

"""
    _harris_priester_mod_density_poly(F10ₐ::Number, i::Integer, c::Integer) -> Number

Evaluate the cubic polynomial in the flux `F10ₐ` [sfu] whose coefficients are stored in
the columns `c` to `c + 3` of the `i`-th row of `_HARRIS_PRIESTER_MOD_COEFS`, returning
the related density [g / km³]. Notice that we build the coefficient tuple explicitly since
splatting a runtime view into `@evalpoly` leads to dynamic dispatch and allocations.
"""
function _harris_priester_mod_density_poly(F10ₐ::Number, i::Integer, c::Integer)
    C = _HARRIS_PRIESTER_MOD_COEFS
    return evalpoly(F10ₐ, (C[i, c], C[i, c + 1], C[i, c + 2], C[i, c + 3]))
end

"""
    _scale_height_junk(x₁::Number, x₂::Number, h::Number, H₁::Number, H₂::Number) -> Number

Blend the scale heights `H₁` and `H₂` [km] defined at the borders `x₁` and `x₂` [km] of the
blending interval at the altitude `h` [km] using the Junkins / Jancaitis weighting method
for third-order continuity [1, Section 3.2].
"""
function _scale_height_junk(x₁::Number, x₂::Number, h::Number, H₁::Number, H₂::Number)
    ξ = (h - x₁) / (x₂ - x₁)

    # Weighting function for third-order continuity [1, Eq. 20].
    w = ξ^4 * @evalpoly(ξ, 35, -84, 70, -20)

    return H₁ + w * (H₂ - H₁)
end
