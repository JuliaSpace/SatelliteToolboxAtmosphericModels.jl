## Description #############################################################################
#
# Harris-Priester atmospheric density model.
#
## References ##############################################################################
#
# [1] Harris, I., Priester, W (1962). Time-dependent structure of the upper atmosphere.
#     Journal of the Atmospheric Sciences, 19(4), pp. 286-301.
#
# [2] Montenbruck, O., Gill, E (2000). Satellite Orbits: Models, Methods and Applications.
#     Springer-Verlag, Berlin, Germany, Section 3.5.1.
#
# [3] Long, A. C., Cappellari Jr., J. O., Velez, C. E., Fuchs, A. J (editors) (1989).
#     Goddard Trajectory Determination System (GTDS) Mathematical Theory (Revision 1).
#     FDD/552-89/0001 and CSC/TR-89/6001, Section 4.4.
#
############################################################################################

export harrispriester

"""
    harrispriester(
        instant::DateTime,
        ϕ_gd::Number,
        λ::Number,
        h::Number;
        kwargs...
    ) -> Number
    harrispriester(jd::Number, ϕ_gd::Number, λ::Number, h::Number; kwargs...) -> Number

Compute the atmospheric density [kg / m³] using the Harris-Priester model [1] with the
diurnal bulge modeling and the density profile for the mean solar activity of [2, 3].

The angle between the diurnal bulge apex and the position is computed using the geodetic
latitude `ϕ_gd` in place of the geocentric declination, as in [2]. The difference is lower
than 0.2° and its effect on the density is negligible for the accuracy of the model.

The model is valid only inside the altitude range of the density profile `alt_ρ` (100 km to
1000 km for the default profile). The function throws an `ArgumentError` if the altitude `h`
is outside this range or if the profile has less than two rows, less than three columns, or
altitudes that are not sorted in ascending order.

# Arguments

- `instant::DateTime`: Instant to compute the model represented using `DateTime`.
- `jd::Number`: Julian day to compute the model.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Geodetic longitude [rad].
- `h::Number`: Geodetic altitude [m].

# Keywords

- `n::Number`: Cosine exponent in the diurnal bulge modeling (`2 <= n <= 7`). If `n` is
    `2`, it models a smooth transition, whereas if `n` is `6`, it models a sharp transition.
    The function throws an `ArgumentError` if `n` is outside this range.
    (**Default**: `4`)
- `alt_ρ::AbstractMatrix`: Matrix containing the minimum and maximum density profiles,
    where the columns are the altitude [m], the minimum density [kg / m³], and the maximum
    density [kg / m³], respectively.
    (**Default**: `_HARRIS_PRIESTER_ALT_RHO`, related to the mean solar activity)

# Returns

- `Number`: Atmospheric density [kg / m³]. The type is the promotion of the types of the
    numeric inputs.

# References

- **[1]** Harris, I., Priester, W (1962). *Time-dependent structure of the upper
    atmosphere*. **Journal of the Atmospheric Sciences**, 19(4), pp. 286-301.
- **[2]** Montenbruck, O., Gill, E (2000). *Satellite Orbits: Models, Methods and
    Applications*. **Springer-Verlag**, Berlin, Germany, Section 3.5.1.
- **[3]** Long, A. C., Cappellari Jr., J. O., Velez, C. E., Fuchs, A. J (editors) (1989).
    *Goddard Trajectory Determination System (GTDS) Mathematical Theory (Revision 1)*.
    FDD/552-89/0001 and CSC/TR-89/6001, Section 4.4.
"""
function harrispriester(
    instant::DateTime,
    ϕ_gd::Number,
    λ::Number,
    h::Number;
    n::Number = 4,
    alt_ρ::AbstractMatrix{<:Number} = _HARRIS_PRIESTER_ALT_RHO,
)
    return harrispriester(datetime2julian(instant), ϕ_gd, λ, h; n = n, alt_ρ = alt_ρ)
end

function harrispriester(
    jd::JT,
    ϕ_gd::PT,
    λ::LT,
    h::HT;
    n::Number = 4,
    alt_ρ::AbstractMatrix{DT} = _HARRIS_PRIESTER_ALT_RHO,
) where {JT <: Number, PT <: Number, LT <: Number, HT <: Number, DT <: Number}
    (2 <= n <= 7) || throw(ArgumentError("The cosine exponent must be between 2 and 7."))

    # Validate the density profile table and the altitude.
    ((size(alt_ρ, 1) >= 2) && (size(alt_ρ, 2) >= 3)) || throw(
        ArgumentError("The density profile must have at least 2 rows and 3 columns.")
    )

    _harris_priester_is_sorted(alt_ρ) || throw(
        ArgumentError("The altitudes in the density profile must be in ascending order.")
    )

    _check_altitude(h, alt_ρ[begin, 1], alt_ρ[end, 1])

    RT = float(promote_type(JT, PT, LT, HT, DT))

    # Compute the Sun declination, the Sun right ascension, and the right ascension of the
    # selected location [rad].
    δs, Ωs, Ωp = _sun_geometry(jd, λ)

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

    # The `max` protects the square root against negative arguments caused by rounding
    # errors when `cos_ψ` is close to -1. The power is skipped on the night side, where it
    # underflows to 0.
    c2ψ2 = max((1 + cos_ψ) / 2, zero(cos_ψ))
    cψ2 = RT(√(c2ψ2))
    cos_pow = (cψ2 > eps(RT)) ? cψ2^n : zero(RT)

    # Search for the altitude index in the density table.
    ia = _harris_priester_layer(alt_ρ, h)

    # Fractional satellite height.
    h₁ = alt_ρ[ia, 1]
    h₂ = alt_ρ[ia + 1, 1]
    dh = (h₁ - h) / (h₁ - h₂)

    # Minimum exponential density interpolation.
    ρ_min₁ = alt_ρ[ia, 2]
    ρ_min₂ = alt_ρ[ia + 1, 2]
    ρ_min = ρ_min₁ * (ρ_min₂ / ρ_min₁)^dh

    iszero(cos_pow) && return RT(ρ_min)

    # Maximum exponential density interpolation.
    ρ_max₁ = alt_ρ[ia, 3]
    ρ_max₂ = alt_ρ[ia + 1, 3]
    ρ_max = ρ_max₁ * (ρ_max₂ / ρ_max₁)^dh
    return RT(ρ_min + (ρ_max - ρ_min) * cos_pow)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _harris_priester_is_sorted(alt_ρ::AbstractMatrix) -> Bool

Return `true` if the altitudes in the first column of the density profile `alt_ρ` are in
ascending order, or `false` otherwise.
"""
function _harris_priester_is_sorted(alt_ρ::AbstractMatrix)
    for k in (firstindex(alt_ρ, 1) + 1):lastindex(alt_ρ, 1)
        (alt_ρ[k, 1] < alt_ρ[k - 1, 1]) && return false
    end

    return true
end

"""
    _harris_priester_layer(alt_ρ::AbstractMatrix, h::Number) -> Int

Return the row index `ia` of the density profile `alt_ρ` such that the altitude `h` [m] is
inside the layer between the rows `ia` and `ia + 1`, using a binary search over the first
column. The altitude must be inside the range of the profile.
"""
function _harris_priester_layer(alt_ρ::AbstractMatrix, h::Number)
    lo = firstindex(alt_ρ, 1)
    hi = lastindex(alt_ρ, 1)

    # Binary search for the last row whose altitude is not higher than `h`.
    while hi - lo > 1
        mid = (lo + hi) >> 1

        if alt_ρ[mid, 1] <= h
            lo = mid
        else
            hi = mid
        end
    end

    return lo
end
