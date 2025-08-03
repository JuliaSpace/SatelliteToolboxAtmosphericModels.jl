module SatelliteToolboxAtmosphericModelsForwardDiffExt

import SatelliteToolboxAtmosphericModels.AtmosphericModels: _jr1971_roots
using ForwardDiff
using ImplicitDifferentiation
using PolynomialRoots

export _jr1971_roots

#==
This file defines the overload of the function `SatelliteToolboxAtmosphericModels.AtmosphericModels._jr1971_roots`
to work with `ForwardDiff.Dual` numbers. This is required to compute the derivatives of the
JR71 model using ForwardDiff.jl.

The problem is that the root finding algorithm (`PolynomialRoots.roots`) is not compatible with
`ForwardDiff.Dual`. To overcome this, we use `ImplicitDifferentiation.jl`.
==#

"""
    _differentiable_roots_forward(x) -> NTuple{2, Vector}

This is the forward pass for the implicit function. It computes the roots from the
polynomial coefficients `x`. The roots are returned as a flat vector of real numbers
(Re, Im, Re, Im, ...).

"""
function _differentiable_roots_forward(x)
    # We need to get the values from the dual numbers to use in `roots`.
    coeffs = ForwardDiff.value.(x)
    r = roots(coeffs; polish = true)

    # We need to sort the roots to have a canonical representation and avoid
    # discontinuities when the order changes. We sort by real part, then
    # imaginary part.
    sort!(r, by = z -> (real(z), imag(z)))

    # Reshape the 4 complex roots into a vector of 8 real numbers.
    y = vcat([[real(root), imag(root)] for root in r]...)
    return y, nothing
end

"""
    _differentiable_roots_conditions(x, y, z) -> Vector

This is the conditions function for the implicit function. It must evaluate to
zero if `y` is the correct set of roots for the polynomial with coefficients
`x`.

"""
function _differentiable_roots_conditions(x, y, z)
    # Reconstruct the complex roots from the real vector `y`.
    roots_complex = [y[1] + im * y[2], y[3] + im * y[4], y[5] + im * y[6], y[7] + im * y[8]]

    # The condition is that P(root) = 0 for each root. This gives us 4 complex
    # equations, which we split into 8 real equations.
    conditions_vec = vcat(
        [[real(@evalpoly(root, x...)), imag(@evalpoly(root, x...))] for root in roots_complex]...
    )
    return conditions_vec
end

# Create the implicit function that computes the roots of a polynomial in a
# differentiable way.
const differentiable_roots = ImplicitFunction(
    _differentiable_roots_forward,
    _differentiable_roots_conditions
)

"""
    _jr1971_roots(p::AbstractVector{<:ForwardDiff.Dual})

Overload of `_jr1971_roots` for `ForwardDiff.Dual` numbers.

This function calls the differentiable root finding algorithm and then
post-processes the results to return the values in the format expected by the
JR71 model.

"""
function _jr1971_roots(p::AbstractVector{<:ForwardDiff.Dual})
    roots_f = differentiable_roots(p)[1]

    # Reconstruct the complex roots with dual numbers.
    r = [roots_f[1] + im * roots_f[2], roots_f[3] + im * roots_f[4], roots_f[5] + im * roots_f[6], roots_f[7] + im * roots_f[8]]

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

    return r₁[1], r₂[1], x[1], y[1]
end

end
