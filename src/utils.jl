## Description #############################################################################
#
# Utility functions shared by the atmospheric models.
#
############################################################################################

"""
    _get_doy(jd::Number) -> Number

Return the fractional day of the year for the Julian day `jd` [UTC], i.e. the number of
days elapsed since the beginning of the year plus one (January 1st at 00:00 is 1.0).
"""
function _get_doy(jd::Number)
    # Compute the year given the selected Julian Day.
    year, _, _, = jd_to_date(jd)

    # Compute the day of the year.
    doy = jd - date_to_jd(year, 1, 1, 0, 0, 0) + 1

    return doy
end

"""
    _sun_geometry(jd::Number, λ::Number) -> Number, Number, Number

Compute the Sun geometry used by the atmospheric models at the Julian day `jd` [UTC] for a
location with longitude `λ` [rad], using the Sun position represented in the MOD reference
frame.

# Returns

- `Number`: Sun declination [rad].
- `Number`: Sun right ascension [rad].
- `Number`: Right ascension of the selected location [rad], i.e. the longitude added to the
    Greenwich mean sidereal time.
"""
function _sun_geometry(jd::Number, λ::Number)
    # Compute the Sun position represented in the inertial reference frame (MOD).
    s_i = sun_position_mod(jd)

    # Compute the Sun declination [rad].
    δs = atan(s_i[3], √(s_i[1] * s_i[1] + s_i[2] * s_i[2]))

    # Compute the Sun right ascension [rad].
    Ωs = atan(s_i[2], s_i[1])

    # Compute the right ascension of the selected location w.r.t. the inertial reference
    # frame.
    Ωp = λ + jd_to_gmst(jd)

    return δs, Ωs, Ωp
end

"""
    _f10_81day_mean(series::Val, jd::Number) -> Float64

Compute the 81-day average of the daily F10.7 solar flux index [sfu] centered on the Julian
day `jd` [UTC], using the space index `series` (`Val(:F10obs)` or `Val(:F10adj)`) fetched
with **SpaceIndices.jl**.
"""
function _f10_81day_mean(series::Val, jd::Number)
    return sum(space_index(series, jd + k) for k in -40:40) / 81
end

"""
    _kp_3h(instant::DateTime) -> Float64

Return the Kp geomagnetic index of the 3-hour interval containing `instant` [UTC], fetched
with **SpaceIndices.jl**. The index is considered constant inside each 3-hour interval.
"""
function _kp_3h(instant::DateTime)
    # Obtain the Kp vector, containing the Kp values for every 3 hours of the day.
    Kp_vect = space_index(Val(:Kp), instant)

    # Get the number of seconds elapsed since the beginning of the day and select the
    # related 3-hour interval.
    day = Date(instant) |> DateTime
    Δt  = Dates.value(instant - day) / 1000
    id  = clamp(floor(Int, Δt / 10_800) + 1, 1, 8)

    return Kp_vect[id]
end
