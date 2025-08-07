## Description #############################################################################
#
# Mathematical functions used in the NRLMSISE-00 model.
#
## References ##############################################################################
#
# [1] https://www.brodo.de/space/nrlmsise/index.html
#
############################################################################################

"""
    _spline_∫(x::NTuple{N, T}, y::NTuple{N, T}, ∂²y::NTuple{N, T}, xf::Number) where {N, T<:Number} -> float(T)

Compute the integral of the cubic spline function `y(x)` from `x[1]` to `xf`, where the
function second derivatives evaluated at `x` are `∂²y`.

# Arguments

- `x::NTuple{N, T}`: X components of the tabulated function in ascending order.
- `y::NTuple{N, T}`: Y components of the tabulated function evaluated at `x`.
- `∂²y::NTuple{N, T}`: Second derivatives of `y(x)` `∂²y/∂x²` evaluated at `x`.
- `xf::Number`: Abscissa endpoint for integration.
"""
function _spline_∫(
    x::NTuple{N, T},
    y::NTuple{N, T},
    ∂²y::NTuple{N, T},
    xf::T
) where {N, T<:Number}

    int = T(0)
    k₀  = 1
    k₁  = 2

    @inbounds while (xf > x[k₀]) && (k₁ <= N)
        xᵢ = ((k₁ <= (N - 1)) && (xf >= x[k₁])) ? x[k₁] : xf
        h  = (x[k₁] - x[k₀])
        a  = (x[k₁] - xᵢ   ) / h
        b  = (xᵢ    - x[k₀]) / h

        a² = a^2
        b² = b^2
        a⁴ = a^4.0
        b⁴ = b^4.0
        h² = h^2

        k_a = -(1 + a⁴) / 4 + a² / 2
        k_b = b⁴ / 4 - b² / 2

        int += h * ((1 - a²) * y[k₀] / 2 + b² *  y[k₁] / 2 + (k_a  * ∂²y[k₀] + k_b * ∂²y[k₁]) * h² / 6)

        k₀ += 1
        k₁ += 1
    end

    return int
end


"""
    _spline_∂²(x::NTuple{N, T}, y::NTuple{N, T}, ∂²y₁::T, ∂²yₙ::T) where {N, T<:Number} -> NTuple{N, T}

Compute the 2nd derivatives of the cubic spline interpolation `y(x)` given the 2nd
derivatives at `x[1]` (`∂²y₁`) and at `x[N]` (`∂²yₙ`). This functions return a tuple with
the evaluated 2nd derivatives at each point in `x`.

!!! note
    This function was adapted from Numerical Recipes.

!!! note
    Values higher than `0.99e30` in the 2nd derivatives at the borders (`∂²y₁` and `∂²yₙ`)
    are interpreted as `0`.

# Arguments

- `x::NTuple{N, T}`: X components of the tabulated function in ascending order.
- `y::NTuple{N, T}`: Y components of the tabulated function evaluated at `x`.
- `∂²y₁::T`: Second derivative of `y(x)` `∂²y/∂x²` evaluated at `x[1]`.
- `∂²yₙ::T`: Second derivative of `y(x)` `∂²y/∂x²` evaluated at `x[N]`.
"""
function _spline_∂²(
    x::NTuple{N, T},
    y::NTuple{N, T},
    ∂²y₁::T,
    ∂²yₙ::T
) where {N, T<:Number}

    # Initialize the tuples that holds the derivatives. Notice that we do not use vectors to
    # avoid allocations. This code can be much slower for very large `N`. However, N <= 5
    # for the NRLMSISE-00 model.
    u   = ntuple(_ -> T(0), Val(N))
    ∂²y = ntuple(_ -> T(0), Val(N))

    if (∂²y₁ > 0.99e30)
        @reset ∂²y[1] = T(0)
        @reset u[1]   = T(0)
    else
        @reset ∂²y[1] = -T(1 / 2)
        @reset u[1]   = (3 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - ∂²y₁)
    end

    @inbounds for i in 2:N-1
        σ             = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        p             = σ * ∂²y[i-1] + 2
        @reset ∂²y[i] = (σ - 1) / p
        m_a           = (y[i+1] - y[i]  ) / (x[i+1] - x[i]  )
        m_b           = (y[i]   - y[i-1]) / (x[i]   - x[i-1])
        @reset u[i]   = (6 * (m_a - m_b) / (x[i+1] - x[i-1]) - σ * u[i-1]) / p
    end

    if ∂²yₙ > 0.99e30
        qₙ = T(0)
        uₙ = T(0)
    else
        qₙ = T(1 / 2)
        uₙ = (3 / (x[N] - x[N-1])) * (∂²yₙ - (y[N] - y[N-1]) / (x[N] - x[N-1]))
    end

    @reset ∂²y[N] = (uₙ - qₙ * u[N-1]) / (qₙ * ∂²y[N-1] + 1)

    @inbounds for k in N-1:-1:1
        @reset ∂²y[k] = ∂²y[k] * ∂²y[k+1] + u[k]
    end

    return ∂²y
end

"""
    _spline(x::NTuple{N, T}, y::NTuple{N, T}, ∂²y::NTuple{N, T}, xᵢ::T) where {N, T<:Number} -> float(T)

Compute the interpolation of the cubic spline `y(x)` with second derivatives `∂²y` at `xᵢ`.

!!! note
    This function was adapted from Numerical Recipes.

# Arguments

- `x::NTuple{N, T}`: X components of the tabulated function in ascending order.
- `y::NTuple{N, T}`: Y components of the tabulated function evaluated at `x`.
- `∂²y::NTuple{N, T}`: Second derivatives of `y(x)` `∂²y/∂x²` evaluated at `x`.
- `xᵢ::T`: Point to compute the interpolation.
"""
function _spline(
    x::NTuple{N, T},
    y::NTuple{N, T},
    ∂²y::NTuple{N, T},
    xᵢ::T
) where {N, T<:Number}

    k₀ = 1
    k₁ = N

    @inbounds while (k₁ - k₀) > 1
        k = div(k₁ + k₀, 2, RoundNearest)

        if x[k] > xᵢ
            k₁ = k
        else
            k₀ = k
        end
    end

    h = x[k₁] - x[k₀]

    (h == 0) &&
        throw(ArgumentError("It is not allowed to have two points with the same abscissa."))

    a  = (x[k₁] - xᵢ   ) / h
    b  = (xᵢ    - x[k₀]) / h
    yᵢ = a * y[k₀] + b * y[k₁] + ((a^3.0 - a) * ∂²y[k₀] + (b^3.0 - b) * ∂²y[k₁]) * h^2 / 6

    return yᵢ
end

"""
    _ζ(r_lat::Number, zz::Number, zl::Number) -> Number

Compute the zeta function.
"""
function _ζ(r_lat::Number, zz::Number, zl::Number)
    return (zz - zl) * (r_lat + zl) / (r_lat + zz)
end
