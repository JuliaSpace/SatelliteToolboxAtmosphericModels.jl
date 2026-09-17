## Description #############################################################################
#
# Types related to the Jacchia 1977 atmospheric model.
#
############################################################################################

export Jacchia1977Output

"""
    struct Jacchia1977Output{T <: Number} <: AbstractAtmosphericModelOutput

Output of the atmospheric model Jacchia 1977.

# Fields

- `total_density::T`: Total atmospheric density [kg / m³].
- `temperature::T`: Local temperature at the selected position [K], computed from the
    temperature profile related to the local quiet exospheric temperature (phase angle of
    -60°, as prescribed for the actual temperature in eq. 26 of the report) increased by
    the geomagnetic variation.
- `exospheric_temperature::T`: Mean exospheric temperature `T½` above the selected position
    [K], as defined in eq. 20 of the report.
- `N2_number_density::T`: Number density of N₂ [1 / m³].
- `O2_number_density::T`: Number density of O₂ [1 / m³].
- `O_number_density::T`: Number density of O [1 / m³].
- `Ar_number_density::T`: Number density of Ar [1 / m³].
- `He_number_density::T`: Number density of He [1 / m³].
- `H_number_density::T`: Number density of H [1 / m³]. It is 0 at altitudes up to 140 km,
    where the model does not include the hydrogen.
"""
struct Jacchia1977Output{T <: Number} <: AbstractAtmosphericModelOutput
    total_density::T
    temperature::T
    exospheric_temperature::T
    N2_number_density::T
    O2_number_density::T
    O_number_density::T
    Ar_number_density::T
    He_number_density::T
    H_number_density::T
end

_model_name(::Type{<:Jacchia1977Output}) = "Jacchia 1977"
_model_description(::Type{<:Jacchia1977Output}) = "Jacchia 1977"

function _show_fields(::Type{<:Jacchia1977Output})
    return (
        ("Total density", :total_density, _SHOW_FORMAT_DENSITY, "kg / m³"),
        ("Temperature", :temperature, _SHOW_FORMAT_TEMPERATURE, "K"),
        ("Exospheric Temp.", :exospheric_temperature, _SHOW_FORMAT_TEMPERATURE, "K"),
        ("N₂ number density", :N2_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("O₂ number density", :O2_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("O  number density", :O_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("Ar number density", :Ar_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("He number density", :He_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
        ("H  number density", :H_number_density, _SHOW_FORMAT_DENSITY, "1 / m³"),
    )
end
