## Description #############################################################################
#
#   Exponential atmosphere model.
#
## Reference ###############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
############################################################################################

export exponential

"""
    exponential(h::Number) -> Float64

Compute the atmospheric density [kg / m³] at the altitude `h` [m] above the ellipsoid using
the exponential atmospheric model:

                    ┌            ┐
                    │    h - h₀  │
    ρ(h) = ρ₀ . exp │ - ──────── │ ,
                    │      H     │
                    └            ┘

in which `ρ₀`, `h₀`, and `H` are parameters obtained from tables that depend only on `h`.
"""
function exponential(h::Number)
    # Check the bounds.
    h < 0 && throw(ArgumentError("The height must be non-negative."))

    # Transform `h` to km.
    h /= 1000

    # Get the values for the exponential model. Since the altitude table is sorted, we can
    # use a binary search to find the layer related to the altitude `h`.
    id  = clamp(searchsortedlast(_EXPONENTIAL_ATMOSPHERE_H₀, h), 1, 28)
    h₀  = _EXPONENTIAL_ATMOSPHERE_H₀[id]
    ρ₀  = _EXPONENTIAL_ATMOSPHERE_ρ₀[id]
    H   = _EXPONENTIAL_ATMOSPHERE_H[id]

    # Compute the density.
    density = ρ₀ * exp(-(h - h₀) / H)

    return density
end
