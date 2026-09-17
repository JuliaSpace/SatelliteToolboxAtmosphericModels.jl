## Description #############################################################################
#
# Exponential atmosphere model.
#
## Reference ###############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
############################################################################################

export exponential

"""
    exponential(h::Number) -> Number

Compute the atmospheric density [kg / m³] at the altitude `h` [m] above the ellipsoid using
the exponential atmospheric model, returning a value with the floating-point type of `h`:

                    ┌            ┐
                    │    h - h₀  │
    ρ(h) = ρ₀ . exp │ - ──────── │ ,
                    │      H     │
                    └            ┘

in which `ρ₀`, `h₀`, and `H` are parameters obtained from tables that depend only on `h`.

The function throws an `ArgumentError` if the altitude `h` is negative. Above 1000 km, the
parameters of the last layer of the table are used.
"""
function exponential(h::Number)
    _check_altitude(h, 0, Inf)

    RT = float(typeof(h))

    # Transform `h` to km.
    h_km = RT(h / 1000)

    # Get the values for the exponential model. Since the altitude table is sorted, we can
    # use a binary search to find the layer related to the altitude `h`. Notice that the
    # index is always valid since `h_km >= 0` and the first altitude of the table is 0.
    id = searchsortedlast(_EXPONENTIAL_ATMOSPHERE_H₀, h_km)
    h₀ = RT(_EXPONENTIAL_ATMOSPHERE_H₀[id])
    ρ₀ = RT(_EXPONENTIAL_ATMOSPHERE_ρ₀[id])
    H  = RT(_EXPONENTIAL_ATMOSPHERE_H[id])

    # Compute the density.
    return ρ₀ * exp(-(h_km - h₀) / H)
end
