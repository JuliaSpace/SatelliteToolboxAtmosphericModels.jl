## Description #############################################################################
#
# Utility functions.
#
############################################################################################

"""
    _get_doy(jd::Number) -> Number

Return the day of the year for the Julian day `jd`.
"""
function _get_doy(jd::Number)
    # Compute the year given the selected Julian Day.
    year, _, _, = jd_to_date(jd)

    # Compute the day of the year.
    doy = jd - date_to_jd(year, 1, 1, 0, 0, 0) + 1

    return doy
end
