## Description #############################################################################
#
# Types related to the Jacchia 1977 atmospheric model.
#
############################################################################################

export Jacchia1977Output

"""
    struct Jacchia1977Output{T<:Number}

Output of the atmospheric model Jacchia 1977.

# Fields

- `total_density::T`: Total atmospheric density [kg / m³].
- `temperature::T`: Temperature at the selected position [K].
- `exospheric_temperature::T`: Mean exospheric temperature `T½` above the selected position
    [K].
- `N2_number_density::T`: Number density of N₂ [1 / m³].
- `O2_number_density::T`: Number density of O₂ [1 / m³].
- `O_number_density::T`: Number density of O [1 / m³].
- `Ar_number_density::T`: Number density of Ar [1 / m³].
- `He_number_density::T`: Number density of He [1 / m³].
- `H_number_density::T`: Number density of H [1 / m³].
"""
struct Jacchia1977Output{T <: Number}
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
