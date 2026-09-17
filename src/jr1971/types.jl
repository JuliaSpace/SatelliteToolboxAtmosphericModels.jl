## Description #############################################################################
#
# Types related to the Jacchia-Roberts 1971 atmospheric model.
#
############################################################################################

export JR1971Output

"""
    struct JR1971Output{T <: Number} <: AbstractAtmosphericModelOutput

Output of the atmospheric model Jacchia-Roberts 1971.

# Fields

- `total_density::T`: Total atmospheric density [kg / m³].
- `temperature::T`: Temperature at the selected position [K].
- `exospheric_temperature::T`: Exospheric temperature [K].
- `N2_number_density::T`: Number density of N₂ [1 / m³].
- `O2_number_density::T`: Number density of O₂ [1 / m³].
- `O_number_density::T`: Number density of O [1 / m³].
- `Ar_number_density::T`: Number density of Ar [1 / m³].
- `He_number_density::T`: Number density of He [1 / m³].
- `H_number_density::T`: Number density of H [1 / m³].
"""
struct JR1971Output{T <: Number} <: AbstractAtmosphericModelOutput
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

_model_name(::Type{<:JR1971Output}) = "JR1971"
_model_description(::Type{<:JR1971Output}) = "Jacchia-Roberts 1971"

function _show_fields(::Type{<:JR1971Output})
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
