## Description #############################################################################
#
# Auxiliary private functions for the NRLMSISE-00 model.
#
############################################################################################

"""
    _ccor(h::T, r::T, h₁::T, zh::T) where T<:Number -> T

Compute the chemistry / dissociation correction for MSIS models.

# Arguments

- `h::Number`: Altitude.
- `r::Number`: Target ratio.
- `h₁::Number`: Transition scale length.
- `zh::Number`: Altitude of `1/2 r`.
"""
function _ccor(h::HT, r::T, h₁::HT2, zh::ZT) where {HT<:Number, T<:Number, HT2<:Number, ZT<:Number}

    RT = promote_type(HT, T, HT2, ZT)

    e = (h - zh) / h₁

    (e > +70) && return exp(RT(0))
    (e < -70) && return exp(r)

    return exp(r / (1 + exp(e)))
end

"""
    _ccor2(alt::T, r::T, h₁::T, zh::T, h₂::T) where T<:Number -> T

Compute the O and O₂ chemistry / dissociation correction for MSIS models.

# Arguments

- `h::Number`: Altitude.
- `r::Number`: Target ration.
- `h₁::Number`: Transition scale length.
- `zh::Number`: Altitude of `1/2 r`.
- `h₂::Number`: Transition scale length 2.
"""
function _ccor2(h::HT, r::T, h₁::HT2, zh::ZT, h₂::HT3) where {HT<:Number, T<:Number, HT2<:Number, ZT<:Number, HT3<:Number}
    e1 = (h - zh) / h₁
    e2 = (h - zh) / h₂

    ((e1 > +70) || (e2 > +70)) && return exp(T(0))
    ((e1 < -70) && (e2 < -70)) && return exp(r)

    return exp(r / (1 + (exp(e1) + exp(e2)) / 2))
end

"""
  _dnet(dd::T, dm::T, zhm::T, xmm::T, xm::T) where T<:Number -> T

Compute the turbopause correction for MSIS models, returning the combined density.

# Arguments

- `dd::T`: Diffusive density.
- `dm::T`: Full mixed density.
- `zhm::T`: Transition scale length.
- `xmm::T`: Full mixed molecular weight.
- `xm::T`: Species molecular weight.
"""
function _dnet(dd::DT, dm::DT2, zhm::ZT, xmm::XT, xm::XT2) where {DT<:Number, DT2<:Number, ZT<:Number, XT<:Number, XT2<:Number}

    RT = promote_type(DT, DT2, ZT, XT, XT2)

    a  = zhm / (xmm - xm)

    if !((dm > 0) && (dd > 0))
        if (dd == 0) && (dm == 0)
            dd = RT(1)
        end

        (dm == 0) && return dd
        (dd == 0) && return dm
    end

    ylog = a * log(dm / dd)

    (ylog < -10) && return dd
    (ylog > +10) && return dm

    return dd * (1 + exp(ylog))^(1 / a)
end

"""
    _gravity_and_effective_radius(ϕ_gd::T) where T<:Number -> T, T

Compute the gravity [cm / s²] and effective radius [km] at the geodetic latitude `ϕ_gd` [°].
"""
function _gravity_and_effective_radius(ϕ_gd::T) where T<:Number
    # Auxiliary variable to reduce computational burden.
    c_2ϕ  = cos(2 * T(_DEG_TO_RAD) * ϕ_gd)

    # Compute the gravity at the selected latitude.
    g_lat = T(_REFERENCE_GRAVITY) * (1 - T(0.0026373) * c_2ϕ)

    # Compute the effective radius at the selected latitude.
    r_lat  = 2 * g_lat / (T(3.085462e-6) + T(2.27e-9) * c_2ϕ) * T(1e-5)

    return g_lat, r_lat
end

"""
    _scale_height(h::T, xm::T, temp::T, g_lat::T, r_lat::T) where T<:Number -> T

Compute the scale height.

# Arguments

- `h::T`: Altitude [km].
- `xm::T`: Species molecular weight [ ].
- `temp::T`: Temperature [K].
- `g_lat::T`: Reference gravity at desired latitude [cm / s²].
- `r_lat::T`: Reference radius at desired latitude [km].
"""
function _scale_height(h::HT, xm::XT, temp::TT, g_lat::GT, r_lat::RLT) where {HT<:Number, XT<:Number, TT<:Number, GT<:Number, RLT<:Number}

    RT = promote_type(HT, XT, TT, GT, RLT)

    # Compute the gravity at the selected altitude.
    gh = g_lat / (1 + h / r_lat)^2

    # Compute the scale height
    sh = RT(_RGAS) * temp / (gh * xm)

    return sh
end

############################################################################################
#                             3hr Magnetic Activity Functions                              #
############################################################################################

"""
    _g0(a::Number, p::AbstractVector)

Compute `g₀` function (see Eq. A24d) using the coefficients `abs_p25 = abs(p[25])` and
`p26 = p[26]`.
"""
function _g₀(a::Number, abs_p25::Number, p26::Number)
    return (a - 4 + (p26 - 1) * (a - 4 + (exp(-abs_p25 * (a - 4)) - 1) / abs_p25))
end

"""
    _sg₀(ex::Number, ap::AbstractVector, abs_p25::Number, p26::Number)

Compute the `sg₀` function (see Eq. A24a) using the `ap` vector and the coefficients
`abs_p25` and `p26`.
"""
function _sg₀(ex::Number, ap::AbstractVector, abs_p25::Number, p26::Number)
    # Auxiliary variables.
    g₂ = _g₀(ap[2], abs_p25, p26)
    g₃ = _g₀(ap[3], abs_p25, p26)
    g₄ = _g₀(ap[4], abs_p25, p26)
    g₅ = _g₀(ap[5], abs_p25, p26)
    g₆ = _g₀(ap[6], abs_p25, p26)
    g₇ = _g₀(ap[7], abs_p25, p26)

    ex²  = ex   * ex
    ex³  = ex²  * ex
    ex⁴  = ex²  * ex²
    ex⁸  = ex⁴  * ex⁴
    ex¹² = ex⁸  * ex⁴
    ex¹⁹ = ex¹² * ex⁴ * ex³

    sumex = 1 + (1 - ex¹⁹) / (1 - ex) * √ex
    r = g₂ + g₃ * ex + g₄ * ex² + g₅ * ex³ + (g₆ * ex⁴ + g₇ * ex¹²) * (1 - ex⁸) / (1 - ex)

    return r / sumex
end
