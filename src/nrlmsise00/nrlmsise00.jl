## Description #############################################################################
#
# The NRLMSISE-00 empirical atmosphere model was developed by Mike Picone, Alan Hedin, and
# Doug Drob based on the MSISE90 model.
#
# The MSISE90 model describes the neutral temperature and densities in Earth's atmosphere
# from ground to thermospheric heights. Below 72.5 km the model is primarily based on the
# MAP Handbook (Labitzke et al., 1985) tabulation of zonal average temperature and pressure
# by Barnett and Corney, which was also used for the CIRA-86. Below 20 km these data were
# supplemented with averages from the National Meteorological Center (NMC). In addition,
# pitot tube, falling sphere, and grenade sounder rocket measurements from 1947 to 1972 were
# taken into consideration. Above 72.5 km MSISE-90 is essentially a revised MSIS-86 model
# taking into account data derived from space shuttle flights and newer incoherent scatter
# results. For someone interested only in the thermosphere (above 120 km), the author
# recommends the MSIS-86 model.  MSISE is also not the model of preference for specialized
# tropospheric work.  It is rather for studies that reach across several atmospheric
# boundaries.
#
#                 (quoted from http://nssdc.gsfc.nasa.gov/space/model/atmos/nrlmsise00.html)
#
# This Julia version of NRLMSISE-00 was converted from the C version implemented and
# maintained by Dominik Brodowski <devel@brodo.de> and available at
# http://www.brodo.de/english/pub/nrlmsise/index.html .
#
# The source code is available at the following git:
#
#     https://git.linta.de/?p=~brodo/nrlmsise-00.git;a=tree
#
# The conversion also used information available at the FORTRAN source code available at
#
#     https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/
#
## References ##############################################################################
#
# [1] https://www.brodo.de/space/nrlmsise/index.html
# [2] https://www.orekit.org/site-orekit-11.0/xref/org/orekit/models/earth/atmosphere/NRLMSISE00.html
#
############################################################################################

export nrlmsise00

"""
    nrlmsise00(instant::DateTime, h::Number, ϕ_gd::Number, λ::Number[, F10ₐ::Number, F10::Number, ap::Union{Number, AbstractVector}]; kwargs...) -> Nrlmsise00Output{Float64}
    nrlmsise00(jd::Number, h::Number, ϕ_gd::Number, λ::Number[, F10ₐ::Number, F10::Number, ap::Union{Number, AbstractVector}]; kwargs...) -> Nrlmsise00Output{Float64}

Compute the atmospheric density using the NRLMSISE-00 model.

If we omit all space indices, the system tries to obtain them automatically for the selected
day `jd` or `instant`. However, the indices must be already initialized using the function
`SpaceIndices.init()`.

# Arguments

- `instant::DateTime`: Instant to compute the model represent using `DateTime`.
- `jd::Number`: Julian day to compute the model.
- `h::Number`: Altitude [m].
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `F10ₐ::Number`: 10.7-cm averaged solar flux, 90-day centered on input time [sfu].
- `F10::Number`: 10.7-cm solar flux [sfu].
- `ap::Union{Number, AbstractVector}`: Magnetic index, see the section **AP** for more
    information.

# Keywords

- `flags::Nrlmsise00Flags`: A list of flags to configure the model. For more information,
    see [`Nrlmsise00Flags`]@(ref). (**Default** = `Nrlmsise00Flags()`)
- `include_anomalous_oxygen::Bool`: If `true`, the anomalous oxygen density will be included
    in the total density computation. (**Default** = `true`)
- `P::Union{Nothing, Matrix}`: If the user passes a matrix with dimensions equal to or
    greater than 8 × 4, it will be used when computing the Legendre associated functions,
    reducing allocations and improving the performance. If it is `nothing`, the matrix is
    allocated inside the function. (**Default** `nothing`)

# Returns

- `Nrlmsise00Output{Float64}`: Structure containing the results obtained from the model.

# AP

The input variable `ap` contains the magnetic index. It can be a `Number` or an
`AbstractVector`.

If `ap` is a number, it must contain the daily magnetic index.

If `ap` is an `AbstractVector`, it must be a vector with 7 dimensions as described below:

| Index | Description                                                                   |
|-------|:------------------------------------------------------------------------------|
|     1 | Daily AP.                                                                     |
|     2 | 3 hour AP index for current time.                                             |
|     3 | 3 hour AP index for 3 hours before current time.                              |
|     4 | 3 hour AP index for 6 hours before current time.                              |
|     5 | 3 hour AP index for 9 hours before current time.                              |
|     6 | Average of eight 3 hour AP indices from 12 to 33 hours prior to current time. |
|     7 | Average of eight 3 hour AP indices from 36 to 57 hours prior to current time. |


# Extended Help

1. The densities of `O`, `H`, and `N` are set to `0` below `72.5 km`.
2. The exospheric temperature is set to global average for altitudes below `120 km`. The
    `120 km` gradient is left at global average value for altitudes below `72.5 km`.
3. Anomalous oxygen is defined as hot atomic oxygen or ionized oxygen that can become
    appreciable at high altitudes (`> 500 km`) for some ranges of inputs, thereby affection
    drag on satellites and debris. We group these species under the term **Anomalous
    Oxygen**, since their individual variations are not presently separable with the drag
    data used to define this model component.

## Notes on Input Variables

`F10` and `F10ₐ` values used to generate the model correspond to the 10.7 cm radio flux at
the actual distance of the Earth from the Sun rather than the radio flux at 1 AU. The
following site provides both classes of values:

    ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/

`F10`, `F10ₐ`, and `ap` effects are neither large nor well established below 80 km and
these parameters should be set to 150, 150, and 4 respectively.

If `include_anomalous_oxygen` is `false`, the `total_density` field in the output is the sum
of the mass densities of the species `He`, `O`, `N₂`, `O₂`, `Ar`, `H`, and `N`, but **does
not** include anomalous oxygen.

If `include_anomalous_oxygen` is `false`, the `total_density` field in the output is the
effective total mass density for drag and is the sum of the mass densities of all species in
this model **including** the anomalous oxygen.
"""
function nrlmsise00(
    instant::DateTime,
    h::Number,
    ϕ_gd::Number,
    λ::Number;
    flags::Nrlmsise00Flags = Nrlmsise00Flags(),
    include_anomalous_oxygen::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    verbose::Val{verbosity} = Val(true),
) where {verbosity}
    return nrlmsise00(
        datetime2julian(instant),
        h,
        ϕ_gd,
        λ;
        flags = flags,
        include_anomalous_oxygen = include_anomalous_oxygen,
        P = P,
        verbose = Val(verbosity),
    )
end

function nrlmsise00(
    jd::JT,
    h::HT,
    ϕ_gd::PT,
    λ::LT;
    flags::Nrlmsise00Flags = Nrlmsise00Flags(),
    include_anomalous_oxygen::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    verbose::Val{verbosity} = Val(true),
) where {JT<:Number, HT<:Number, PT<:Number, LT<:Number, verbosity}

    # Fetch the space indices.
    #
    # NOTE: If the altitude is lower than 80km, set to default according to the instructions
    # in NRLMSISE-00 source code.
    if h < 80e3
        F10ₐ = 150.0
        F10  = 150.0
        ap   = 4.0

        verbosity && @debug """
        NRLMSISE00 - Using default indices since h < 80 km
          Daily F10.7           : $(F10) sfu
          90-day avareged F10.7 : $(F10ₐ) sfu
          Ap                    : $(ap)
        """
    else
        # TODO: The online version of NRLMSISE-00 seems to use 89 days, whereas the
        # NRLMSISE-00 source code mentions 81 days.
        F10ₐ = sum((space_index.(Val(:F10adj), jd + k) for k in -45:44)) / 90
        F10  = space_index(Val(:F10adj), jd - 1)
        ap   = sum(space_index(Val(:Ap), jd)) / 8

        verbosity && @debug """
        NRLMSISE00 - Fetched Space Indices
          Daily F10.7           : $(F10) sfu
          90-day avareged F10.7 : $(F10ₐ) sfu
          Ap                    : $(ap)
        """
    end

    # Call the NRLMSISE-00 model.
    return nrlmsise00(
        jd,
        h,
        ϕ_gd,
        λ,
        F10ₐ,
        F10,
        ap;
        flags = flags,
        include_anomalous_oxygen = include_anomalous_oxygen,
        P = P
    )
end

function nrlmsise00(
    instant::DateTime,
    h::Number,
    ϕ_gd::Number,
    λ::Number,
    F10ₐ::Number,
    F10::Number,
    ap::Union{Number, AbstractVector};
    flags::Nrlmsise00Flags = Nrlmsise00Flags(),
    include_anomalous_oxygen::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing
)
    return nrlmsise00(
        datetime2julian(instant),
        h,
        ϕ_gd,
        λ,
        F10ₐ,
        F10,
        ap;
        flags = flags,
        include_anomalous_oxygen = include_anomalous_oxygen,
        P = P
    )
end

function nrlmsise00(
    jd::JT,
    h::HT,
    ϕ_gd::PT,
    λ::LT,
    F10ₐ::FT,
    F10::FT2,
    ap::T_AP;
    flags::Nrlmsise00Flags = Nrlmsise00Flags(),
    include_anomalous_oxygen::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing
) where {JT<:Number, HT<:Number, PT<:Number, LT<:Number, FT<:Number, FT2<:Number, T_AP<:Union{Number, AbstractVector}}

    RT = promote_type(JT, HT, PT, LT, FT, FT2)

    # == Compute Auxiliary Variables =======================================================

    # Convert the Julian Day to Date.
    Y, _, _, _, _, _ = jd_to_date(jd)

    # Get the number of days since the beginning of the year.
    doy = _get_doy(jd)

    # Get the number of seconds since the beginning of the day.
    Δds = jd - floor(jd)

    # Get the local apparent solar time [hours].
    #
    # TODO: To be very precise, I think this should also take into consideration the
    # Equation of Time. However, the online version of NRLMSISE-00 does not use this.
    lst = Δds / 3600 + λ * 12 / π

    df  = F10  - F10ₐ
    dfa = F10ₐ - 150

    stloc,  ctloc  = sincos(1 * _HOUR_TO_RAD * lst)
    s2tloc, c2tloc = sincos(2 * _HOUR_TO_RAD * lst)
    s3tloc, c3tloc = sincos(3 * _HOUR_TO_RAD * lst)

    # Compute Legendre polynomials.
    #
    # Notice that the original NRLMSISE-00 algorithms considers that the Legendre matrix is
    # upper triangular, whereas we use the lower triangular representation. Hence, we need
    # to transpose it.
    #
    # Furthermore, the NRLMSISE-00 algorithm only uses terms with maximum degree 7 and
    # maximum order 3.
    if isnothing(P)
        plg = legendre(Val(:unnormalized), π / 2 - ϕ_gd, 7, 3; ph_term = false)'
    else
        rows, cols = size(P)

        if (rows < 8) || (cols < 4)
            throw(ArgumentError("The matrix P must have at least 8 × 4 elements."))
        end

        legendre!(Val(:unnormalized), P, π / 2 - ϕ_gd, 7, 3; ph_term = false)
        plg = P'
    end

    # == Latitude Variation of Gravity =====================================================
    #
    # None for flags.time_independent = false.
    g_lat, r_lat = _gravity_and_effective_radius(
        (!flags.time_independent) ? RT(_REFERENCE_LATITUDE) : RT(ϕ_gd / _DEG_TO_RAD)
    )

    # == Create the NRLMSISE00 Structure ===================================================

    nrlmsise00d = Nrlmsise00Structure{RT, T_AP}(
        Y,
        doy,
        Δds,
        h / 1000,
        ϕ_gd / _DEG_TO_RAD,
        λ / _DEG_TO_RAD,
        lst,
        F10ₐ,
        F10,
        ap,
        flags,
        r_lat,
        g_lat,
        df,
        dfa,
        plg,
        ctloc,
        stloc,
        c2tloc,
        s2tloc,
        c3tloc,
        s3tloc,
        0,
        0,
        0,
        0,
        0
    )

    # Call the NRLMSISE-00 model.
    ~, nrlmsise00_out = include_anomalous_oxygen ?
        _gtd7d(nrlmsise00d) :
        _gtd7(nrlmsise00d)

    return nrlmsise00_out
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _densm(h::T, d0::T, xm::T, tz::T, r_lat::T, g_lat::T, tn2::NTuple{N2, T}, tgn2::NTuple{2, T}, tn3::NTuple{N3, T}, tgn3::NTuple{2, T}) where {N2<:Interger, N3<:Integer, T<:Number} -> float(T), float(T)

Compute the temperature and density profiles for the lower atmosphere.

!!! note
    This function returns the density if `xm` is not 0, or the temperature otherwise.

# Arguments

- `h::T`: Altitude [km].
- `d₀::T`: Reference density, returned if `h > _ZN2[1]`.
- `xm::T`: Species molecular weight [ ].
- `g_lat::T`: Reference gravity at desired latitude [cm / s²].
- `r_lat::T`: Reference radius at desired latitude [km].
- `tn2::NTuple{N2, T}`: Temperature at the nodes for ZN2 scale [K].
- `tgn2::NTuple{N2, T}`: Temperature gradients at the end nodes for ZN2 scale.
- `tn3::NTuple{N3, T}`: Temperature at the nodes for ZN3 scale [K].
- `tgn3::NTuple{N3, T}`: Temperature gradients at the end nodes for ZN3 scale.

# Returns

- `T`: Density [1 / cm³] is `xm` is not 0, or the temperature [K] otherwise.
"""
function _densm(
    h::HT,
    d₀::DT,
    xm::XT,
    g_lat::GT,
    r_lat::RLT,
    tn2::NTuple{4, TNT},
    tgn2::NTuple{2, TGT},
    tn3::NTuple{5, TNT2},
    tgn3::NTuple{2, TGT2},
) where {HT<:Number, DT<:Number, XT<:Number, GT<:Number, RLT<:Number, TNT<:Number, TGT<:Number, TNT2<:Number, TGT2<:Number}

    RT = promote_type(HT, DT, XT, GT, RLT, TNT, TGT, TNT2, TGT2)

    # == Initialization of Variables =======================================================

    density = d₀
    (h > _ZN2[begin]) && return density

    # == Stratosphere / Mesosphere Temperature =============================================

    z     = (h > _ZN2[end]) ? h : _ZN2[end]
    z1    = _ZN2[begin]
    z2    = _ZN2[end]
    t1    = tn2[begin]
    t2    = tn2[end]
    zg    = _ζ(r_lat, z, z1)
    zgdif = _ζ(r_lat, z2, z1)

    # Set up spline nodes.
    xs2 = ntuple(_ -> RT(0), Val(4))
    ys2 = ntuple(_ -> RT(0), Val(4))

    @inbounds for k in 1:4
        @reset xs2[k] = _ζ(r_lat, _ZN2[k], z1) / zgdif
        @reset ys2[k] = 1 / tn2[k]
    end

    ∂²y₁ = -tgn2[begin] / (t1 * t1) * zgdif
    ∂²yₙ = -tgn2[end]   / (t2 * t2) * zgdif * ((r_lat + z2) / (r_lat + z1))^2

    # Calculate spline coefficients.
    ∂²y = _spline_∂²(xs2, ys2, ∂²y₁, ∂²yₙ)

    # Interpolate at desired point.
    x = zg / zgdif
    y = _spline(xs2, ys2, ∂²y, x)

    # Temperature at altitude.
    tz = 1 / y

    if xm != 0
        # Compute the gravity at `z1`.
        g_h = g_lat / (1 + z1 / r_lat)^2

        # Calculate stratosphere / mesosphere density.
        γ = xm * g_h * zgdif / RT(_RGAS)

        # Integrate temperature profile.
        expl = min(γ * _spline_∫(xs2, ys2, ∂²y, x), RT(50));

        # Density at altitude.
        density *= (t1 / tz) * exp(-expl)
    end

    if h > _ZN3[1]
        if xm == 0
            return tz
        else
            return density
        end
    end

    # == Troposhepre / Stratosphere Temperature ============================================

    z     = h
    z1    = RT(_ZN3[begin])
    z2    = RT(_ZN3[end])
    t1    = tn3[begin]
    t2    = tn3[end]
    zg    = _ζ(r_lat, z, z1)
    zgdif = _ζ(r_lat, z2, z1)

    # Set up spline nodes.
    xs3 = ntuple(_ -> RT(0), Val(5))
    ys3 = ntuple(_ -> RT(0), Val(5))

    @inbounds for k in 1:5
        @reset xs3[k] = _ζ(r_lat, _ZN3[k], z1) / zgdif
        @reset ys3[k] = 1 / tn3[k]
    end

    ∂²y₁ = -tgn3[begin] / (t1 * t1) * zgdif
    ∂²yₙ = -tgn3[end]   / (t2 * t2) * zgdif * ((r_lat + z2) / (r_lat + z1))^2

    # Calculate spline coefficients.
    ∂²y = _spline_∂²(xs3, ys3, ∂²y₁, ∂²yₙ)
    x   = zg / zgdif
    y   = _spline(xs3, ys3, ∂²y, x)

    # Temperature at altitude.
    tz = 1 / y

    if xm != 0
        # Compute the gravity at `z1`.
        g_h = g_lat / (1 + z1 / r_lat)^2

        # Calculate tropospheric / stratosphere density.
        γ = xm * g_h * zgdif / RT(_RGAS);

        # Integrate temperature profile.
        expl = min(γ * _spline_∫(xs3, ys3, ∂²y, x) , RT(50))

        # Density at altitude.
        density *= (t1 / tz) * exp(-expl)

        return density
    else
        return tz
    end
end

"""
    _densu(h::T, dlb::T, tinf::T, tlb::T, xm::T, α::T, zlb::T, s2::T, g_lat::T, r_lat::T, tn1::NTuple{5, T}, tgn1::NTuple{2, T}) where T<:Number -> T, NTuple{5, T}, NTuple{2, T}

Compute the density [1 / cm³] or temperature [K] profiles according to the new lower thermo
polynomial.

!!! note
    This function returns the density if `xm` is not 0, or the temperature otherwise.

# Arguments

- `h::T`: Altitude [km].
- `dlb::T`: Density at lower boundary [1 / cm³].
- `tinf::T`: Exospheric temperature [K].
- `tlb::T`: Temperature at lower boundary [K].
- `xm::T`: Species molecular weight [ ].
- `α::T`: Thermal diffusion coefficient.
- `zlb::T`: Altitude at lower boundary [km].
- `s2::T`: Slope.
- `g_lat::T`: Reference gravity at the latitude [cm / s²].
- `r_lat::T`: Reference radius at the latitude [km].
- `tn1::NTuple{5, T}`: Temperature at nodes for ZN1 scale [K].
- `tgn1::NTuple{2, T}`: Temperature gradients at end nodes for ZN1 scale.

# Returns

- `T`: Density [1 / cm³] is `xm` is not 0, or the temperature [K] otherwise.
- `NTuple{5, T}`: Updated `tn1`.
- `NTuple{2, T}`: Updated `tgn1`.
"""
function _densu(
    h::HT,
    dlb::DT,
    tinf::TinfT,
    tlb::TT,
    xm::XT,
    α::AT,
    zlb::ZT,
    s2::ST,
    g_lat::GLT,
    r_lat::RLT,
    tn1::NTuple{5, TNT},
    tgn1::NTuple{2, TGT}
) where {HT<:Number, DT<:Number, TinfT<:Number, TT<:Number, XT<:Number, AT<:Number, ZT<:Number, ST<:Number, GLT<:Number, RLT<:Number, TNT<:Number, TGT<:Number}

    RT = promote_type(HT, DT, TinfT, TT, XT, AT, ZT, ST, GLT, RLT, TNT, TGT)

    x     = RT(0)
    z1    = RT(0)
    t1    = RT(0)
    zgdif = RT(0)
    xs    = ntuple(_ -> RT(0), Val(5))
    ys    = ntuple(_ -> RT(0), Val(5))
    ∂²y   = ntuple(_ -> RT(0), Val(5))

    # Joining altitudes of Bates and spline.
    z = max(h, _ZN1[begin])

    # Geopotential altitude difference from ZLB.
    zg2 = _ζ(r_lat, z, zlb)

    # Bates temperature.
    tt = tinf - (tinf - tlb) * exp(-s2 * zg2)
    ta = tt
    tz = tt

    @inbounds if h < _ZN1[begin]
        # Compute the temperature below ZA temperature gradient at ZA from Bates profile.
        dta = (tinf - ta) * s2 * ((r_lat + zlb) / (r_lat + _ZN1[begin]))^2

        @reset tgn1[begin] = RT(dta)
        @reset tn1[begin]  = RT(ta)
        z  = (h > _ZN1[end]) ? h : _ZN1[end]
        z1 = _ZN1[begin]
        z2 = _ZN1[end]
        t1 = tn1[begin]
        t2 = tn1[end]

        # Geopotential difference from z1.
        zg    = _ζ(r_lat,  z, z1)
        zgdif = _ζ(r_lat, z2, z1)

        # Set up spline nodes.
        for k in 1:5
            @reset xs[k] = _ζ(r_lat, _ZN1[k], z1) / zgdif
            @reset ys[k] = 1 / tn1[k]
        end

        # End node derivatives.
        ∂²y₁ = -tgn1[begin] / (t1 * t1) * zgdif
        ∂²yₙ = -tgn1[end]   / (t2 * t2) * zgdif * ((r_lat + z2) / (r_lat + z1))^2

        # Compute spline coefficients.
        @reset ∂²y = _spline_∂²(xs, ys, ∂²y₁, ∂²yₙ)

        # Interpolate at the desired point.
        x = zg / zgdif
        y = _spline(xs, ys, ∂²y, x)

        # Temperature at altitude.
        tz = 1 / y
    end

    (xm == 0) && return tz, tn1, tgn1

    # Compute the gravity at `zlb`.
    g_h = g_lat / (1 + zlb / r_lat)^2

    # Calculate density above _ZN1[1].
    γ = xm * g_h / (s2 * RT(_RGAS) * tinf)
    expl = exp(-s2 * γ * zg2)

    if (expl > 50) || (tt <= 0)
        expl = RT(50)
    end

    # Density at altitude.
    density = dlb * (tlb / tt)^(1 + α + γ) * expl

    (h >= _ZN1[1]) && return density, tn1, tgn1

    # Compute the gravity at `z1`.
    g_h = g_lat / (1 + z1 / r_lat)^2

    # Compute density below _ZN1[1].
    γ = xm * g_h * zgdif / RT(_RGAS)

    # Integrate spline temperatures.
    expl = γ * _spline_∫(xs, ys, ∂²y, x)

    if (expl > 50) || (tz <= 0)
        expl = RT(50)
    end

    # Density at altitude.
    density *= (t1 / tz)^(1 + α) * exp(-expl)

    return density, tn1, tgn1
end

"""
    _globe7(nrlmsise00d::Nrlmsise00Structure{T}, p::AbstractVector{T}) where T<:Number -> Nrlmsise00Structure{T}, T

Compute the function `G(L)` with upper thermosphere parameters `p` and the NRLMSISE-00
structure `nrlmsise00`.

!!! note
    The variables `apt` and `apdf` inside `nrlmsise00d` can be modified inside this
    function.

# Returns

- `Nrlmsise00Structure{T}`: Modified structure `nrlmsise00d`.
- `T`: Result of `G(L)`.
"""
function _globe7(nrlmsise00d::Nrlmsise00Structure{T}, p::AbstractVector{V}) where {T<:Number, V<:Number}
    # == Unpack NRLMSISE00 Structure =======================================================

    ap     = nrlmsise00d.ap
    apdf   = nrlmsise00d.apdf
    apt    = nrlmsise00d.apt
    c2tloc = nrlmsise00d.c2tloc
    c3tloc = nrlmsise00d.c3tloc
    ctloc  = nrlmsise00d.ctloc
    df     = nrlmsise00d.df
    dfa    = nrlmsise00d.dfa
    doy    = nrlmsise00d.doy
    flags  = nrlmsise00d.flags
    lst    = nrlmsise00d.lst
    plg    = nrlmsise00d.plg
    s2tloc = nrlmsise00d.s2tloc
    s3tloc = nrlmsise00d.s3tloc
    sec    = nrlmsise00d.sec
    stloc  = nrlmsise00d.stloc
    λ      = nrlmsise00d.λ
    ϕ_gd   = nrlmsise00d.ϕ_gd

    # == Initialization of Variables =======================================================

    t₁  = T(0)
    t₂  = T(0)
    t₃  = T(0)
    t₄  = T(0)
    t₅  = T(0)
    t₆  = T(0)
    t₇  = T(0)
    t₈  = T(0)
    t₉  = T(0)
    t₁₀ = T(0)
    t₁₁ = T(0)
    t₁₂ = T(0)
    t₁₃ = T(0)
    t₁₄ = T(0)

    tloc = lst

    cd32 = cos(1 * _DAY_TO_RAD * (doy - p[32]))
    cd18 = cos(2 * _DAY_TO_RAD * (doy - p[18]))
    cd14 = cos(1 * _DAY_TO_RAD * (doy - p[14]))
    cd39 = cos(2 * _DAY_TO_RAD * (doy - p[39]))

    # == F10.7 Effect ======================================================================

    t₁ = p[20] * df *(1 + p[60] * dfa) + p[21] * df^2 + p[22] * dfa + p[30] * dfa^2
    f1 = 1 + (p[48] * dfa + p[20] * df + p[21] * df^2) * flags.F10_Mean
    f2 = 1 + (p[50] * dfa + p[20] * df + p[21] * df^2) * flags.F10_Mean

    # == Time Independent ==================================================================

    t₂ = p[2]  * plg[1, 3] + p[3] * plg[1, 5] + p[23] * plg[1, 7] + p[27] * plg[1, 2] +
         p[15] * plg[1, 3] * dfa * flags.F10_Mean

    # == Symmetrical Annual ================================================================

    t₃ = p[19] * cd32

    # == Symmetrical Semiannual ============================================================

    t₄ = (p[16] + p[17] * plg[1, 3]) * cd18

    # == Asymmetrical Annual ===============================================================

    t₅ = f1 * (p[10] * plg[1, 2] + p[11] * plg[1, 4]) * cd14

    # == Asymmetrical Semiannual ===========================================================

    t₆ = p[38] * plg[1, 2] * cd39

    # == Diurnal ===========================================================================

    if flags.diurnal
        t71 = (p[12] * plg[2, 3]) * cd14 * flags.asym_annual
        t72 = (p[13] * plg[2, 3]) * cd14 * flags.asym_annual

        t₇ = f2 * (
            (p[4] * plg[2, 2] + p[5] * plg[2, 4] + p[28] * plg[2, 6] + t71) * ctloc +
            (p[7] * plg[2, 2] + p[8] * plg[2, 4] + p[29] * plg[2, 6] + t72) * stloc
        )
    end

    # == Semidiurnal =======================================================================

    if flags.semidiurnal
        t81 = (p[24] * plg[3, 4] + p[36] * plg[3, 6]) * cd14 * flags.asym_annual
        t82 = (p[34] * plg[3, 4] + p[37] * plg[3, 6]) * cd14 * flags.asym_annual

        t₈ = f2 * (
            (p[6] * plg[3, 3] + p[42] * plg[3, 5] + t81) * c2tloc +
            (p[9] * plg[3, 3] + p[43] * plg[3, 5] + t82) * s2tloc
        )
    end

    # == Terdiurnal ========================================================================

    if flags.terdiurnal
        t91 = (p[94] * plg[4, 5] + p[47] * plg[4, 7]) * cd14 * flags.asym_annual
        t92 = (p[95] * plg[4, 5] + p[49] * plg[4, 7]) * cd14 * flags.asym_annual

        t₁₄ = f2 * ((p[40] * plg[4, 4] + t91) * s3tloc + (p[41] * plg[4, 4] + t92) * c3tloc)
    end

    # == Magnetic Activity Based on Daily AP ===============================================

    if ap isa AbstractVector
        if p[52] != 0
            exp1 = min(exp(-10800 * abs(p[52]) / (1 + p[139] * (45 - abs(ϕ_gd)))), 0.99999)

            if p[25] < 1.0e-4
                p[25] = 1.0e-4
            end

            apt = _sg₀(exp1, ap, abs(p[25]), p[26])
            aux = cos(_HOUR_TO_RAD * (tloc - p[132]))

            t₉ = apt * (
                (p[51] + p[97] * plg[1, 3] + p[55] * plg[1, 5]) +
                (p[126] * plg[1, 2] + p[127] * plg[1, 4] + p[128] * plg[1, 6]) * cd14 * flags.asym_annual +
                (p[129] * plg[2, 2] + p[130] * plg[2, 4] + p[131] * plg[2, 6]) * aux * flags.diurnal
            )
        end

    else
        apd = ap - 4
        p44 = p[44]
        p45 = p[45]

        if p44 < 0
            p44 = 1e-5
        end

        apdf = apd + (p45 - 1) * (apd + (exp(-p44 * apd) - 1) / p44)

        if flags.daily_ap
            aux = cos(_HOUR_TO_RAD * (tloc - p[125]))

            t₉ = apdf * (
                (p[33] + p[46] * plg[1, 3] + p[35] * plg[1, 5]) +
                (p[101] * plg[1, 2] + p[102] * plg[1, 4] + p[103] * plg[1, 6]) * cd14 * flags.asym_annual +
                (p[122] * plg[2, 2] + p[123] * plg[2, 4] + p[124] * plg[2, 6]) * aux * flags.diurnal
            )
        end
    end

    if flags.all_ut_long_effects && (λ > - 1000)

        # == Longitudinal ==================================================================

        if flags.longitudinal
            sin_g_long, cos_g_long = sincos(_DEG_TO_RAD * λ)

            k₁ = p[65]  * plg[2, 3] + p[66]  * plg[2, 5] + p[67]  * plg[2, 7] +
                 p[104] * plg[2, 2] + p[105] * plg[2, 4] + p[106] * plg[2, 6]
            k₂ = p[110] * plg[2, 2] + p[111] * plg[2, 4] + p[112] * plg[2, 6]
            k₃ = p[91]  * plg[2, 3] + p[92]  * plg[2, 5] + p[93]  * plg[2, 7] +
                 p[107] * plg[2, 2] + p[108] * plg[2, 4] + p[109] * plg[2, 6]
            k₄ = p[113] * plg[2, 2] + p[114] * plg[2, 4] + p[115] * plg[2, 6]

            t₁₁ = (1 + p[81] * dfa * flags.F10_Mean) * (
                (k₁ + flags.asym_annual * k₂ * cd14) * cos_g_long +
                (k₃ + flags.asym_annual * k₄ * cd14) * sin_g_long
            )
        end

        # == UT and Mixed UT, Longitude ====================================================

        if flags.ut_mixed_ut_long

            k₁ = (1 +  p[96] * plg[1, 2]) * (1 + p[82] * dfa * flags.F10_Mean)
            k₂ = 1 + p[120] * plg[1, 2] * flags.asym_annual * cd14
            k₃ = p[69] * plg[1, 2] + p[70] * plg[1, 4] + p[71] * plg[1, 6]
            k₄ = 1 + p[138] * dfa * flags.F10_Mean
            k₅ = p[77] * plg[3, 4] + p[78] * plg[3, 6] + p[79] * plg[3, 8]

            aux₁ = cos(_SEC_TO_RAD * (sec - p[72]))
            aux₂ = cos(_SEC_TO_RAD * (sec - p[80]) + 2 * _DEG_TO_RAD * λ)

            t₁₂ =  k₁ * k₂ * k₃ * aux₁ + flags.longitudinal * k₄ * k₅ * aux₂
        end

        # == UT, Longitude Magnetic Activity ===============================================

        if flags.mixed_ap_ut_long
            if ap isa AbstractVector
                if p[52] != 0

                    k₁ = p[53]  * plg[2, 3] + p[99]  * plg[2, 5] + p[68]  * plg[2, 7]
                    k₂ = p[134] * plg[2, 2] + p[135] * plg[2, 4] + p[136] * plg[2, 6]
                    k₃ = p[56]  * plg[1, 2] + p[57]  * plg[1, 4] + p[58]  * plg[1, 6]

                    aux₁ = cos(_DEG_TO_RAD * (λ - p[98]))
                    aux₂ = cos(_DEG_TO_RAD * (λ - p[137]))
                    aux₃ = cos(_SEC_TO_RAD * (sec - p[59]))

                    t₁₃ = apt * flags.longitudinal * (1 + p[133] * plg[1, 2]) * k₁ * aux₁ +
                          apt * flags.longitudinal * flags.asym_annual * k₂ * cd14 * aux₂ +
                          apt * flags.ut_mixed_ut_long * k₃ * aux₃
                end
            else
                k₁ = p[61]  * plg[2, 3] + p[62]  * plg[2, 5] + p[63]  * plg[2, 7]
                k₂ = p[116] * plg[2, 2] + p[117] * plg[2, 4] + p[118] * plg[2, 6]
                k₃ = p[84]  * plg[1, 2] + p[85]  * plg[1, 4] + p[86]  * plg[1, 6]

                aux₁ = cos(_DEG_TO_RAD * (λ - p[64]))
                aux₂ = cos(_DEG_TO_RAD * (λ - p[119]))
                aux₃ = cos(_SEC_TO_RAD * (sec - p[76]))

                t₁₃ = apdf * flags.longitudinal * (1 + p[121] * plg[1,2]) * k₁ * aux₁ +
                      apdf * flags.longitudinal * flags.asym_annual * k₂ * cd14 * aux₂ +
                      apdf * flags.ut_mixed_ut_long * k₃ * aux₃
            end
        end
    end

    # Update the NRLMSISE-00 structure.
    @reset nrlmsise00d.apt  = T(apt)
    @reset nrlmsise00d.apdf = T(apdf)

    # Parameters not used: 82, 89, 99, 139-149.
    tinf = p[31] +
           flags.F10_Mean            * t₁  +
           flags.time_independent    * t₂  +
           flags.sym_annual          * t₃  +
           flags.sym_semiannual      * t₄  +
           flags.asym_annual         * t₅  +
           flags.asym_semiannual     * t₆  +
           flags.diurnal             * t₇  +
           flags.semidiurnal         * t₈  +
           flags.daily_ap            * t₉  +
           flags.all_ut_long_effects * t₁₀ +
           flags.longitudinal        * t₁₁ +
           flags.ut_mixed_ut_long    * t₁₂ +
           flags.mixed_ap_ut_long    * t₁₃ +
           flags.terdiurnal          * t₁₄

    return nrlmsise00d, tinf
end

"""
    _glob7s(nrlmsise00d::Nrlmsise00Structure{T}, p::AbstractVector{T}) where T<:Number -> T

Compute the function `G(L)` with lower atmosphere parameters `p` and the NRLMSISE-00
structure `nrlmsise00d`.
"""
function _glob7s(nrlmsise00d::Nrlmsise00Structure{T}, p::AbstractVector{V}) where {T<:Number, V<:Number}

    # == Unpack NRLMSISE00 Structure =======================================================

    ap     = nrlmsise00d.ap
    apdf   = nrlmsise00d.apdf
    apt    = nrlmsise00d.apt
    c2tloc = nrlmsise00d.c2tloc
    c3tloc = nrlmsise00d.c3tloc
    ctloc  = nrlmsise00d.ctloc
    dfa    = nrlmsise00d.dfa
    doy    = nrlmsise00d.doy
    flags  = nrlmsise00d.flags
    plg    = nrlmsise00d.plg
    s2tloc = nrlmsise00d.s2tloc
    s3tloc = nrlmsise00d.s3tloc
    stloc  = nrlmsise00d.stloc
    λ      = nrlmsise00d.λ

    # == Initialization of Variables =======================================================

    t₁  = T(0)
    t₂  = T(0)
    t₃  = T(0)
    t₄  = T(0)
    t₅  = T(0)
    t₆  = T(0)
    t₇  = T(0)
    t₈  = T(0)
    t₉  = T(0)
    t₁₀ = T(0)
    t₁₁ = T(0)
    t₁₂ = T(0)
    t₁₃ = T(0)
    t₁₄ = T(0)

    # Confirm parameter set.
    if p[100] == 0
        p[100] = T(2)
    end

    (p[100] != 2) && error("Wrong parameter set for `_glob7s!`.")

    cd32 = cos(1 * _DAY_TO_RAD * (doy - p[32]))
    cd18 = cos(2 * _DAY_TO_RAD * (doy - p[18]))
    cd14 = cos(1 * _DAY_TO_RAD * (doy - p[14]))
    cd39 = cos(2 * _DAY_TO_RAD * (doy - p[39]))

    # == F10.7 =============================================================================

    t₁ = p[22] * dfa

    # == Time Independent ==================================================================

    t₂ = p[2]  * plg[1, 3] + p[3]  * plg[1, 5] + p[23] * plg[1, 7] + p[27] * plg[1, 2] +
         p[15] * plg[1, 4] + p[60] * plg[1, 6]

    # == Symmetrical Annual ================================================================

    t₃ = (p[19] + p[48] * plg[1, 3] + p[30] * plg[1, 5]) * cd32

    # == Symmetrical Semiannual ============================================================

    t₄ = (p[16] + p[17] * plg[1, 3] + p[31] * plg[1, 5]) * cd18

    # == Asymmetrical Annual ===============================================================

    t₅ = (p[10] * plg[1, 2] + p[11] * plg[1, 4] + p[21] * plg[1, 6]) * cd14

    # == Asymmetrical Semiannual ===========================================================

    t₆ = p[38] * plg[1, 2] * cd39

    # == Diurnal ===========================================================================

    if flags.diurnal
        t71 = p[12] * plg[2, 3] * cd14 * flags.asym_annual
        t72 = p[13] * plg[2, 3] * cd14 * flags.asym_annual
        t₇  = (p[4] * plg[2, 2] + p[5] * plg[2, 4] + t71) * ctloc +
              (p[7] * plg[2, 2] + p[8] * plg[2, 4] + t72) * stloc
    end

    # == Semidiurnal =======================================================================

    if flags.semidiurnal
        t81 = (p[24] * plg[3, 4] + p[36] * plg[3, 6]) * cd14 * flags.asym_annual
        t82 = (p[34] * plg[3, 4] + p[37] * plg[3, 6]) * cd14 * flags.asym_annual
        t₈  = (p[6] * plg[3, 3] + p[42] * plg[3, 5] + t81) * c2tloc +
              (p[9] * plg[3, 3] + p[43] * plg[3, 5] + t82) * s2tloc
    end

    # == Terdiurnal ========================================================================

    if flags.terdiurnal
        t₁₄ = p[40] * plg[4, 4] * s3tloc + p[41] * plg[4, 4] * c3tloc
    end

    # == Magnetic Activity =================================================================

    if flags.daily_ap
            t₉ = p[51] * apt + p[97] * plg[1, 3] * apt * flags.time_independent
        if ap isa AbstractVector
        else
            t₉ = apdf * (p[33] + p[46] * plg[1, 3] * flags.time_independent)
        end
    end

    # == Longitudinal ======================================================================

    if !(!flags.all_ut_long_effects || !flags.longitudinal || (λ <= -1000.0))
        sin_g_long, cos_g_long = sincos(_DEG_TO_RAD * λ)

        k₁ = p[65] * plg[2, 3] + p[66] * plg[2, 5] + p[67] * plg[2, 7] + p[75] * plg[2, 2] +
             p[76] * plg[2, 4] + p[77] * plg[2, 6]
        k₂ = p[91] * plg[2, 3] + p[92] * plg[2, 5] + p[93] * plg[2, 7] + p[78] * plg[2, 2] +
             p[79] * plg[2, 4] + p[80] * plg[2, 6]

        t₁₁ = (k₁ * cos_g_long + k₂ * sin_g_long) * (
            1 + plg[1,2] * (
                p[81] * cos(1 * _DAY_TO_RAD * (doy - p[82])) * flags.asym_annual +
                p[86] * cos(2 * _DAY_TO_RAD * (doy - p[87])) * flags.asym_semiannual
            ) + p[84] * cos(1 * _DAY_TO_RAD * (doy - p[85])) * flags.sym_annual +
                p[88] * cos(2 * _DAY_TO_RAD * (doy - p[89])) * flags.sym_semiannual
        )
    end

    tinf = flags.F10_Mean            * t₁  +
           flags.time_independent    * t₂  +
           flags.sym_annual          * t₃  +
           flags.sym_semiannual      * t₄  +
           flags.asym_annual         * t₅  +
           flags.asym_semiannual     * t₆  +
           flags.diurnal             * t₇  +
           flags.semidiurnal         * t₈  +
           flags.daily_ap            * t₉  +
           flags.all_ut_long_effects * t₁₀ +
           flags.longitudinal        * t₁₁ +
           flags.ut_mixed_ut_long    * t₁₂ +
           flags.mixed_ap_ut_long    * t₁₃ +
           flags.terdiurnal          * t₁₄

    return tinf
end

"""
    _gtd7(nrlmsise00d::Nrlmsise00Structure{T}) where T<:Number -> Nrlmsise00Structure{T}, Nrlmsise00Output{T}

Compute the temperatures and densities using the information inside the structure
`nrlmsise00d` without including the anomalous oxygen in the total density.

# Returns

- `Nrlmsise00Structure{T}`: Modified structure `nrlmsise00d`.
- `Nrlmsise00Output{T}`: Structure with the output information.
"""
function _gtd7(nrlmsise00d::Nrlmsise00Structure{T}) where T<:Number

    # == Constants =========================================================================

    pdm_1  = _NRLMSISE00_PDM_1
    pdm_3  = _NRLMSISE00_PDM_3
    pdm_4  = _NRLMSISE00_PDM_4
    pdm_5  = _NRLMSISE00_PDM_5
    pma_1  = _NRLMSISE00_PMA_1
    pma_2  = _NRLMSISE00_PMA_2
    pma_3  = _NRLMSISE00_PMA_3
    pma_4  = _NRLMSISE00_PMA_4
    pma_5  = _NRLMSISE00_PMA_5
    pma_6  = _NRLMSISE00_PMA_6
    pma_7  = _NRLMSISE00_PMA_7
    pma_8  = _NRLMSISE00_PMA_8
    pma_10 = _NRLMSISE00_PMA_10
    pavgm  = _NRLMSISE00_PAVGM

    # == Unpack NRLMSISE00 Structure =======================================================

    flags = nrlmsise00d.flags
    g_lat = nrlmsise00d.g_lat
    h     = nrlmsise00d.h
    r_lat = nrlmsise00d.r_lat

    # == Initialization of Variables =======================================================

    meso_tn2  = ntuple(_ -> T(0), 4)
    meso_tn3  = ntuple(_ -> T(0), 5)
    meso_tgn2 = ntuple(_ -> T(0), 2)
    meso_tgn3 = ntuple(_ -> T(0), 2)

    # == Latitude Variation of Gravity =====================================================

    xmm = pdm_3[5]

    # == Thermosphere / Mesosphere (above _ZN2[1]) =========================================

    if h < _ZN2[begin]
        @reset nrlmsise00d.h = T(_ZN2[begin])
    end

    nrlmsise00d, out_thermo = _gts7(nrlmsise00d)
    @reset nrlmsise00d.h = h

    # Unpack the values again because `gts7` may have modified `nrlmsise00d`.
    meso_tn1_5  = nrlmsise00d.meso_tn1_5
    meso_tgn1_2 = nrlmsise00d.meso_tgn1_2
    dm28        = nrlmsise00d.dm28

    # If we are above `_ZN2[1]`, then we do not need to compute anything else.
    h >= _ZN2[begin] && return nrlmsise00d, out_thermo

    # Unpack the output values from the thermospheric portion.
    total_density          = out_thermo.total_density
    temperature            = out_thermo.temperature
    exospheric_temperature = out_thermo.exospheric_temperature
    N_number_density       = out_thermo.N_number_density
    N2_number_density      = out_thermo.N2_number_density
    O_number_density       = out_thermo.O_number_density
    aO_number_density      = out_thermo.aO_number_density
    O2_number_density      = out_thermo.O2_number_density
    H_number_density       = out_thermo.H_number_density
    He_number_density      = out_thermo.He_number_density
    Ar_number_density      = out_thermo.Ar_number_density

    # Convert the unit to SI.
    dm28m = dm28 * T(1e6)

    # == Lower Mesosphere / Upper Stratosphere (between `_ZN3[1]` and `_ZN2[1]`) ===========

    @reset meso_tgn2[1] = meso_tgn1_2
    @reset meso_tn2[1]  = meso_tn1_5
    @reset meso_tn2[2]  = pma_1[1] * pavgm[1] / (1 - flags.all_tn2_var * _glob7s(nrlmsise00d, pma_1))
    @reset meso_tn2[3]  = pma_2[1] * pavgm[2] / (1 - flags.all_tn2_var * _glob7s(nrlmsise00d, pma_2))
    @reset meso_tn2[4]  = pma_3[1] * pavgm[3] / (1 - flags.all_tn2_var * flags.all_tn3_var * _glob7s(nrlmsise00d, pma_3))
    @reset meso_tn3[1]  = meso_tn2[4]

    @reset meso_tgn2[2] = pavgm[9] * pma_10[1] * (
        1 + flags.all_tn2_var * flags.all_tn3_var * _glob7s(nrlmsise00d, pma_10)
    ) * meso_tn2[4]^2 / (pma_3[1] * pavgm[3])^2

    # == Lower Stratosphere and Troposphere (below `zn3[1]`) ===============================

    if h < _ZN3[begin]
        @reset meso_tgn3[1] = meso_tgn2[2]
        @reset meso_tn3[2]  = pma_4[1] * pavgm[4] / (1 - flags.all_tn3_var * _glob7s(nrlmsise00d, pma_4))
        @reset meso_tn3[3]  = pma_5[1] * pavgm[5] / (1 - flags.all_tn3_var * _glob7s(nrlmsise00d, pma_5))
        @reset meso_tn3[4]  = pma_6[1] * pavgm[6] / (1 - flags.all_tn3_var * _glob7s(nrlmsise00d, pma_6))
        @reset meso_tn3[5]  = pma_7[1] * pavgm[7] / (1 - flags.all_tn3_var * _glob7s(nrlmsise00d, pma_7))
        @reset meso_tgn3[2] = pma_8[1] * pavgm[8] * (
            1 + flags.all_tn3_var * _glob7s(nrlmsise00d, pma_8)
        ) * meso_tn3[5] * meso_tn3[5] / (pma_7[1] * pavgm[7])^2
    end

    # == Linear Transition to Full Mixing Below `_ZN2[1]` ==================================

    dmc  = (h > _ZMIX) ? 1 - (T(_ZN2[begin]) - h) / (T(_ZN2[begin]) - T(_ZMIX)) : T(0)
    dz28 = N2_number_density

    # == N₂ Density ========================================================================

    dmr = N2_number_density / dm28m - 1

    N2_number_density = _densm(
        h,
        dm28m,
        xmm,
        g_lat,
        r_lat,
        meso_tn2,
        meso_tgn2,
        meso_tn3,
        meso_tgn3
    )

    N2_number_density *= 1 + dmr * dmc

    # == He Density ========================================================================

    dmr = He_number_density / (dz28 * pdm_1[2]) - 1
    He_number_density = N2_number_density * pdm_1[2] * (1 + dmr * dmc)

    # == O Density =========================================================================

    O_number_density = T(0)
    aO_number_density = T(0)

    # == O₂ Density ========================================================================

    dmr = O2_number_density / (dz28 * pdm_4[2]) - 1
    O2_number_density = N2_number_density * pdm_4[2] * (1 + dmr * dmc)

    # == Ar Density ========================================================================

    dmr = Ar_number_density / (dz28 * pdm_5[2]) - 1
    Ar_number_density = N2_number_density * pdm_5[2] * (1 + dmr * dmc)

    # == H Density =========================================================================

    H_number_density = T(0)

    # == N Density =========================================================================

    N_number_density = T(0)

    # == Total Mass Density ================================================================

    total_density = 1.66e-24 * (
        4  * He_number_density +
        16 * O_number_density  +
        28 * N2_number_density +
        32 * O2_number_density +
        40 * Ar_number_density +
        1  * H_number_density  +
        14 * N_number_density
    )

    # Convert the units to SI.
    total_density /= 1000

    # == Temperature at Selected Altitude ==================================================

    temperature = _densm(
        h,
        T(1),
        T(0),
        g_lat,
        r_lat,
        meso_tn2,
        meso_tgn2,
        meso_tn3,
        meso_tgn3
    )

    # Create output structure and return.
    nrlmsise00_out = Nrlmsise00Output{T}(
        total_density,
        temperature,
        exospheric_temperature,
        N_number_density,
        N2_number_density,
        O_number_density,
        aO_number_density,
        O2_number_density,
        H_number_density,
        He_number_density,
        Ar_number_density
    )

    return nrlmsise00d, nrlmsise00_out
end

"""
    _gtd7d(nrlmsise00d::Nrlmsise00Structure{T}) where T<:Number -> Nrlmsise00Structure{T}, Nrlmsise00Output{T}

Compute the temperatures and densities using the information inside the structure
`nrlmsise00d` including the anomalous oxygen in the total density.

# Returns

- `Nrlmsise00Structure{T}`: Modified structure `nrlmsise00d`.
- `Nrlmsise00Output{T}`: Structure with the output information.
"""
function _gtd7d(nrlmsise00d::Nrlmsise00Structure{T}) where T<:Number
    # Call `_gt7d!` to compute the NRLMSISE-00 outputs.
    nrlmsise00d, out = _gtd7(nrlmsise00d)

    # Update the computation of the total mass density.
    total_density = 1.66e-24 * (
        4  * out.He_number_density +
        16 * out.O_number_density  +
        28 * out.N2_number_density +
        32 * out.O2_number_density +
        40 * out.Ar_number_density +
        1  * out.H_number_density  +
        14 * out.N_number_density  +
        16 * out.aO_number_density
    )

    # Convert the unit to SI.
    total_density /= 1000

    # Create the new output and return.
    nrlmsise00_out = Nrlmsise00Output{T}(
        total_density,
        out.temperature,
        out.exospheric_temperature,
        out.N_number_density,
        out.N2_number_density,
        out.O_number_density,
        out.aO_number_density,
        out.O2_number_density,
        out.H_number_density,
        out.He_number_density,
        out.Ar_number_density
    )

    return nrlmsise00d, nrlmsise00_out
end

"""
    _gts7(nrlmsise00d::Nrlmsise00Structure{T}) where T<:Number -> Nrlmsise00Structure{T}, Nrlmsise00Output{T}

Compute the temperatures and densities using the information inside the structure
`nrlmsise00d` and including the anomalous oxygen in the total density for altitudes higher
than 72.5 km (thermospheric portion of NRLMSISE-00).

# Returns

- `Nrlmsise00Structure{T}`: Modified structure `nrlmsise00d`.
- `Nrlmsise00Output{T}`: Structure with the output information.
"""
function _gts7(nrlmsise00d::Nrlmsise00Structure{T}) where T<:Number

    # == Constants =========================================================================

    pdl_1   = _NRLMSISE00_PDL_1
    pdl_2   = _NRLMSISE00_PDL_2
    ptm     = _NRLMSISE00_PTM
    pdm_1   = _NRLMSISE00_PDM_1
    pdm_2   = _NRLMSISE00_PDM_2
    pdm_3   = _NRLMSISE00_PDM_3
    pdm_4   = _NRLMSISE00_PDM_4
    pdm_5   = _NRLMSISE00_PDM_5
    pdm_6   = _NRLMSISE00_PDM_6
    pdm_7   = _NRLMSISE00_PDM_7
    pdm_8   = _NRLMSISE00_PDM_8
    pt      = _NRLMSISE00_PT
    ps      = _NRLMSISE00_PS
    pd_TLB  = _NRLMSISE00_PD_TLB
    pd_N2   = _NRLMSISE00_PD_N2
    ptl_1   = _NRLMSISE00_PTL_1
    ptl_2   = _NRLMSISE00_PTL_2
    ptl_3   = _NRLMSISE00_PTL_3
    ptl_4   = _NRLMSISE00_PTL_4
    pma_9   = _NRLMSISE00_PMA_9
    pd_He   = _NRLMSISE00_PD_HE
    pd_O    = _NRLMSISE00_PD_O
    pd_O2   = _NRLMSISE00_PD_O2
    pd_Ar   = _NRLMSISE00_PD_AR
    pd_H    = _NRLMSISE00_PD_H
    pd_N    = _NRLMSISE00_PD_N
    pd_hotO = _NRLMSISE00_PD_HOT_O

    # Thermal diffusion coefficients for the species.
    α = (-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0)

    # Net density computation altitude limits for the species.
    altl = (200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0)

    # == Unpack NRLMSISE00 structure =======================================================

    dfa   = nrlmsise00d.dfa
    dm28  = nrlmsise00d.dm28
    doy   = nrlmsise00d.doy
    flags = nrlmsise00d.flags
    g_lat = nrlmsise00d.g_lat
    h     = nrlmsise00d.h
    r_lat = nrlmsise00d.r_lat
    ϕ_gd  = nrlmsise00d.ϕ_gd

    # == Initialization of Variables =======================================================

    temperature = T(0)
    meso_tn1  = ntuple(_ -> T(0), 5)
    meso_tgn1 = ntuple(_ -> T(0), 2)

    # == Tinf variations not important below `za` or `zn1[1]` ==============================

    if h > _ZN1[1]
        nrlmsise00d, G_L = _globe7(nrlmsise00d, pt)
        tinf = ptm[1] * pt[1] * (1 + flags.all_tinf_var * G_L)
    else
        tinf = ptm[1] * pt[1]
    end

    exospheric_temperature = tinf

    # == Gradient variations not important below `zn1[5]` ==================================

    if h > _ZN1[5]
        nrlmsise00d, G_L = _globe7(nrlmsise00d, ps)
        g0 = ptm[4] * ps[1] * (1 + flags.all_s_var * G_L)
    else
        g0 = ptm[4] * ps[1]
    end

    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_TLB)
    tlb = ptm[2] * (1 + flags.all_tlb_var * G_L) * pd_TLB[1]
    s   = g0 / (tinf - tlb)

    # Lower thermosphere temperature variations not significant for density above 300 km.

    if h < 300
        @reset meso_tn1[2]  = T(ptm[7] * ptl_1[1] / (1 - flags.all_tn1_var * _glob7s(nrlmsise00d, ptl_1)))
        @reset meso_tn1[3]  = T(ptm[3] * ptl_2[1] / (1 - flags.all_tn1_var * _glob7s(nrlmsise00d, ptl_2)))
        @reset meso_tn1[4]  = T(ptm[8] * ptl_3[1] / (1 - flags.all_tn1_var * _glob7s(nrlmsise00d, ptl_3)))
        @reset meso_tn1[5]  = T(ptm[5] * ptl_4[1] / (1 - flags.all_tn1_var * flags.all_tn2_var * _glob7s(nrlmsise00d, ptl_4)))
        @reset meso_tgn1[2] = T(ptm[9] * pma_9[1] * (
            1 + flags.all_tn1_var * flags.all_tn2_var * _glob7s(nrlmsise00d, pma_9)
        ) * meso_tn1[5]^2 / (ptm[5] * ptl_4[1])^2)
    else
        @reset meso_tn1[2]  = T(ptm[7] * ptl_1[1])
        @reset meso_tn1[3]  = T(ptm[3] * ptl_2[1])
        @reset meso_tn1[4]  = T(ptm[8] * ptl_3[1])
        @reset meso_tn1[5]  = T(ptm[5] * ptl_4[1])
        @reset meso_tgn1[2] = T(ptm[9] * pma_9[1] * meso_tn1[5]^2 / (ptm[5] * ptl_4[1])^2)
    end

    # N2 variation factor at Zlb.
    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_N2)
    g28 = flags.all_nlb_var * G_L

    # == Variation of Turbopause Height ====================================================

    zhf = pdl_2[25] * (
        1 + flags.asym_annual * pdl_1[25] * sin(_DEG_TO_RAD * ϕ_gd) * cos(_DAY_TO_RAD * (doy - pt[14]))
    )
    xmm = pdm_3[5]

    # == N₂ Density ========================================================================

    # Diffusive density at Zlb.
    db28 = pdm_3[1] * exp(g28) * pd_N2[1]

    # Diffusive density at desired altitude.
    N2_number_density, meso_tn1, meso_tgn1 = _densu(
        h,
        db28,
        tinf,
        tlb,
        T(28),
        α[3],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    # Turbopause.
    zh28  = pdm_3[3] * zhf
    zhm28 = pdm_3[4] * pdl_2[6]
    xmd   = 28 - xmm

    # Mixed density at Zlb.
    b28, meso_tn1, meso_tgn1 = _densu(
        zh28,
        db28,
        tinf,
        tlb,
        xmd,
        α[3] - 1,
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    if flags.departures_from_eq && (h < altl[3])
        # Mixed density at desired altitude.
        dm28, meso_tn1, meso_tgn1 = _densu(
            h,
            b28,
            tinf,
            tlb,
            xmm,
            α[3],
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        # Net density at desired altitude.
        N2_number_density = _dnet(N2_number_density, dm28, zhm28, xmm, T(28))
    end

    # == He Density ========================================================================

    # Density variation factor at Zlb.
    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_He)
    g4 = flags.all_nlb_var * G_L

    # Diffusive density at Zlb.
    db04 = pdm_1[1] * exp(g4) * pd_He[1]

    # Diffusive density at desired altitude.
    He_number_density, meso_tn1, meso_tgn1 = _densu(
        h,
        db04,
        tinf,
        tlb,
        T(4),
        α[1],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    if flags.departures_from_eq && (h < altl[3])
        # Turbopause.
        zh04 = pdm_1[3]

        # Mixed density at Zlb.
        b04, meso_tn1, meso_tgn1 = _densu(
            zh04,
            db04,
            tinf,
            tlb,
            4 - xmm,
            α[1] - 1,
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        # Mixed density at desired altitude.
        dm04, meso_tn1, meso_tgn1 = _densu(
            h,
            b04,
            tinf,
            tlb,
            xmm,
            T(0),
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        zhm04 = zhm28

        # Net density at desired altitude.
        He_number_density = _dnet(He_number_density, dm04, zhm04, xmm, T(4))

        # Correction to specified mixing ration at ground.
        rl   = log(b28 * pdm_1[2] / b04)
        zc04 = pdm_1[5] * pdl_2[1]
        hc04 = pdm_1[6] * pdl_2[2]

        # Net density corrected at desired altitude.
        He_number_density *= _ccor(h, rl, hc04, zc04)
    end

    # == O Density =========================================================================

    # Density variation factor at Zlb.
    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_O)
    g16 = flags.all_nlb_var * G_L

    #  Diffusive density at Zlb.
    db16 = pdm_2[1] * exp(g16) * pd_O[1]

    # Diffusive density at desired altitude.
    O_number_density, meso_tn1, meso_tgn1 = _densu(
        h,
        db16,
        tinf,
        tlb,
        T(16),
        α[2],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    if flags.departures_from_eq && (h <= altl[2])
        # Turbopause.
        zh16 = pdm_2[3]

        # Mixed density at Zlb.
        b16, meso_tn1, meso_tgn1 = _densu(
            zh16,
            db16,
            tinf,
            tlb,
            16 - xmm,
            α[2] - 1,
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        # Mixed density at desired altitude.
        dm16, meso_tn1, meso_tgn1 = _densu(
            h,
            b16,
            tinf,
            tlb,
            xmm,
            T(0),
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        zhm16 = zhm28

        # Net density at desired altitude.
        O_number_density  = _dnet(O_number_density, dm16, zhm16, xmm, T(16))
        rl     = pdm_2[2] * pdl_2[17] * (1 + flags.F10_Mean * pdl_1[24] * dfa)
        hc16   = pdm_2[6] * pdl_2[4]
        zc16   = pdm_2[5] * pdl_2[3]
        hc216  = pdm_2[6] * pdl_2[5]
        O_number_density *= _ccor2(h, rl, hc16, zc16, hc216)

        # Chemistry correction.
        hcc16 = pdm_2[8] * pdl_2[14]
        zcc16 = pdm_2[7] * pdl_2[13]
        rc16  = pdm_2[4] * pdl_2[15]

        # Net density corrected at desired altitude.
        O_number_density *= _ccor(h, rc16, hcc16, zcc16)
    end

    # == O₂ Density ========================================================================

    # Density variation factor at Zlb.
    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_O2)
    g32 = flags.all_nlb_var * G_L

    # Diffusive density at Zlb.
    db32 = pdm_4[1] * exp(g32) * pd_O2[1]

    # Diffusive density at desired altitude.
    O2_number_density, meso_tn1, meso_tgn1 = _densu(
        h,
        db32,
        tinf,
        tlb,
        T(32),
        α[4],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    if flags.departures_from_eq
        if h <= altl[4]
            # Turbopause.
            zh32 = pdm_4[3]

            # Mixed density at Zlb.
            b32, meso_tn1, meso_tgn1 = _densu(
                zh32,
                db32,
                tinf,
                tlb,
                32 - xmm,
                α[4] - 1,
                ptm[6],
                s,
                g_lat,
                r_lat,
                meso_tn1,
                meso_tgn1
            )

            # Mixed density at desired altitude.
            dm32, meso_tn1, meso_tgn1 = _densu(
                h,
                b32,
                tinf,
                tlb,
                xmm,
                T(0),
                ptm[6],
                s,
                g_lat,
                r_lat,
                meso_tn1,
                meso_tgn1
            )

            zhm32 = zhm28

            # Net density at desired altitude.
            O2_number_density = _dnet(O2_number_density, dm32, zhm32, xmm, T(32))

            # Correction to specified mixing ratio at ground.
            rl      = log(b28 * pdm_4[2] / b32)
            hc32    = pdm_4[6] * pdl_2[8]
            zc32    = pdm_4[5] * pdl_2[7]
            O2_number_density *= _ccor(h, rl, hc32, zc32)
        end

        # Correction for general departure from diffusive equilibrium above Zlb.
        hcc32  = pdm_4[8] * pdl_2[23]
        hcc232 = pdm_4[8] * pdl_1[23]
        zcc32  = pdm_4[7] * pdl_2[22]
        rc32   = pdm_4[4] * pdl_2[24] * (1 + flags.F10_Mean * pdl_1[24] * dfa)

        # Net density corrected at desired altitude.
        O2_number_density *= _ccor2(h, rc32, hcc32, zcc32, hcc232)
    end

    # == Ar Density ========================================================================

    # Density variation factor at Zlb.
    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_Ar)
    g40 = flags.all_nlb_var * G_L

    # Diffusive density at Zlb.
    db40 = pdm_5[1] * exp(g40) * pd_Ar[1]

    # Diffusive density at desired altitude.
    Ar_number_density, meso_tn1, meso_tgn1= _densu(
        h,
        db40,
        tinf,
        tlb,
        T(40),
        α[5],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    if flags.departures_from_eq && (h <= altl[5])
        # Turbopause.
        zh40 = pdm_5[3]

        # Mixed density at Zlb.
        b40, meso_tn1, meso_tgn1 = _densu(
            zh40,
            db40,
            tinf,
            tlb,
            40 - xmm,
            α[5] - 1,
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        # Mixed density at desired altitude.
        dm40, meso_tn1, meso_tgn1 = _densu(
            h,
            b40,
            tinf,
            tlb,
            xmm,
            T(0),
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        zhm40 = zhm28

        # Net density at desired altitude.
        Ar_number_density = _dnet(Ar_number_density, dm40, zhm40, xmm, T(40))

        # Correction to specified mixing ratio at ground.
        rl   = log(b28 * pdm_5[2] / b40)
        hc40 = pdm_5[6] * pdl_2[10]
        zc40 = pdm_5[5] * pdl_2[9]

        # Net density corrected at desired altitude.
        Ar_number_density *= _ccor(h, rl, hc40, zc40)
    end

    # == H Density =========================================================================

    # Density variation factor at Zlb.
    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_H)
    g1 = flags.all_nlb_var * G_L

    # Diffusive density at Zlb.
    db01 = pdm_6[1] * exp(g1) * pd_H[1]

    # Diffusive density at desired altitude.
    H_number_density, meso_tn1, meso_tgn1 = _densu(
        h,
        db01,
        tinf,
        tlb,
        T(1),
        α[7],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    if flags.departures_from_eq && (h <= altl[7])
        # Turbopause.
        zh01 = pdm_6[3]

        # Mixed density at Zlb.
        b01, meso_tn1, meso_tgn1 = _densu(
            zh01,
            db01,
            tinf,
            tlb,
            1 - xmm,
            α[7] - 1,
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        # Mixed density at desired altitude.
        dm01, meso_tn1, meso_tgn1 = _densu(
            h,
            b01,
            tinf,
            tlb,
            xmm,
            T(0),
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        zhm01 = zhm28

        # Net density at desired altitude.
        H_number_density = _dnet(H_number_density, dm01, zhm01, xmm, T(1))

        # Correction to specified mixing ratio at ground.
        rl     = log(b28 * pdm_6[2] * abs(pdl_2[18]) / b01)
        hc01   = pdm_6[6] * pdl_2[12]
        zc01   = pdm_6[5] * pdl_2[11]
        H_number_density *= _ccor(h, rl, hc01, zc01)

        # Chemistry correction.
        hcc01 = pdm_6[8] * pdl_2[20]
        zcc01 = pdm_6[7] * pdl_2[19]
        rc01  = pdm_6[4] * pdl_2[21]

        # Net density corrected at desired altitude.
        H_number_density *= _ccor(h, rc01, hcc01, zcc01)
    end

    # == N Density =========================================================================

    # Density variation factor at Zlb.
    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_N)
    g14 = flags.all_nlb_var * G_L

    # Diffusive density at Zlb.
    db14 = pdm_7[1] * exp(g14) * pd_N[1]

    # Diffusive density at desired altitude.
    N_number_density, meso_tn1, meso_tgn1 = _densu(
        h,
        db14,
        tinf,
        tlb,
        T(14),
        α[8],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    if flags.departures_from_eq && (h <= altl[8])
        # Turbopause.
        zh14 = pdm_7[3]

        # Mixed density at Zlb.
        b14, meso_tn1, meso_tgn1 = _densu(
            zh14,
            db14,
            tinf,
            tlb,
            14 - xmm,
            α[8] - 1,
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        #  Mixed density at desired altitude.
        dm14, meso_tn1, meso_tgn1 = _densu(
            h,
            b14,
            tinf,
            tlb,
            xmm,
            T(0),
            ptm[6],
            s,
            g_lat,
            r_lat,
            meso_tn1,
            meso_tgn1
        )

        zhm14 = zhm28

        # Net density at desired altitude.
        N_number_density = _dnet(N_number_density, dm14, zhm14, xmm, T(14))

        # Correction to specified mixing ratio at ground.
        rl     = log(b28 * pdm_7[2] * abs(pdl_1[3]) / b14)
        hc14   = pdm_7[6] * pdl_1[2]
        zc14   = pdm_7[5] * pdl_1[1]
        N_number_density *= _ccor(h, rl, hc14, zc14)

        # Chemistry correction.
        hcc14 = pdm_7[8] * pdl_1[5]
        zcc14 = pdm_7[7] * pdl_1[4]
        rc14  = pdm_7[4] * pdl_1[6]

        # Net density corrected at desired altitude.
        N_number_density *= _ccor(h, rc14, hcc14, zcc14)
    end

    # == Anomalous O Density ===============================================================

    nrlmsise00d, G_L = _globe7(nrlmsise00d, pd_hotO)
    g16h  = flags.all_nlb_var * G_L
    db16h = pdm_8[1] * exp(g16h) * pd_hotO[1]
    tho   = pdm_8[10] * pdl_1[7]

    aO_number_density, meso_tn1, meso_tgn1 = _densu(
        h,
        db16h,
        tho,
        tho,
        T(16),
        α[9],
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    zsht = pdm_8[6]
    zmho = pdm_8[5]
    zsho = _scale_height(zmho, T(16), tho, g_lat, r_lat)

    aO_number_density *= exp(-zsht / zsho * (exp(-(h - zmho) / zsht) - 1))

    # == Total Mass Density ================================================================

    total_density = 1.66e-24 * (
        4  * He_number_density +
        16 * O_number_density  +
        28 * N2_number_density +
        32 * O2_number_density +
        40 * Ar_number_density +
        1  * H_number_density  +
        14 * N_number_density
    )

    # == Temperature at Selected Altitude ==================================================

    temperature, meso_tn1, meso_tgn1 = _densu(
        abs(h),
        T(1),
        tinf,
        tlb,
        T(0),
        T(0),
        ptm[6],
        s,
        g_lat,
        r_lat,
        meso_tn1,
        meso_tgn1
    )

    # == Output ============================================================================

    # Convert the result to SI.
    total_density     *= T(1e3)
    N_number_density  *= T(1e6)
    N2_number_density *= T(1e6)
    O_number_density  *= T(1e6)
    aO_number_density *= T(1e6)
    O2_number_density *= T(1e6)
    H_number_density  *= T(1e6)
    He_number_density *= T(1e6)
    Ar_number_density *= T(1e6)

    # Repack variables that were modified.
    @reset nrlmsise00d.meso_tn1_5  = T(meso_tn1[5])
    @reset nrlmsise00d.meso_tgn1_2 = T(meso_tgn1[2])
    @reset nrlmsise00d.dm28        = T(dm28)

    # Create output structure and return.
    nrlmsise00_out = Nrlmsise00Output{T}(
        total_density,
        temperature,
        exospheric_temperature,
        N_number_density,
        N2_number_density,
        O_number_density,
        aO_number_density,
        O2_number_density,
        H_number_density,
        He_number_density,
        Ar_number_density
    )

    return nrlmsise00d, nrlmsise00_out
end
