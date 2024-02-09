## Description #############################################################################
#
# Function to show the results of the NRLMSISE-00 atmospheric model.
#
############################################################################################

function show(io::IO, out::Nrlmsise00Output)
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? string(_B) : ""
    d = color ? string(_D) : ""

    print(io, "$(b)NRLMSISE-00 output$(d) (ρ = ", @sprintf("%g", out.total_density), " kg / m³)")
    return nothing
end

function show(io::IO, mime::MIME"text/plain", out::Nrlmsise00Output)
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? string(_B) : ""
    d = color ? string(_D) : ""

    str_total_density   = @sprintf("%15g", out.total_density)
    str_temperature     = @sprintf("%15.2f", out.temperature)
    str_exospheric_temp = @sprintf("%15.2f", out.exospheric_temperature)
    str_N_number_den    = @sprintf("%15g", out.N_number_density)
    str_N₂_number_den   = @sprintf("%15g", out.N2_number_density)
    str_O_number_den    = @sprintf("%15g", out.O_number_density)
    str_aO_number_den   = @sprintf("%15g", out.aO_number_density)
    str_O₂_number_den   = @sprintf("%15g", out.O2_number_density)
    str_H_number_den    = @sprintf("%15g", out.H_number_density)
    str_He_number_den   = @sprintf("%15g", out.He_number_density)
    str_Ar_number_den   = @sprintf("%15g", out.Ar_number_density)

    println(io, "NRLMSISE-00 Atmospheric Model Result:")
    println(io, "$(b)          Total density :$(d)", str_total_density,   "  kg / m³")
    println(io, "$(b)            Temperature :$(d)", str_temperature,     "  K")
    println(io, "$(b)       Exospheric Temp. :$(d)", str_exospheric_temp, "  K")
    println(io, "$(b)      N  number density :$(d)", str_N_number_den,    "  1 / m³")
    println(io, "$(b)      N₂ number density :$(d)", str_N₂_number_den,   "  1 / m³")
    println(io, "$(b)      O  number density :$(d)", str_O_number_den,    "  1 / m³")
    println(io, "$(b)  Anomalous O num. den. :$(d)", str_aO_number_den,    "  1 / m³")
    println(io, "$(b)      O₂ number density :$(d)", str_O₂_number_den,   "  1 / m³")
    println(io, "$(b)      Ar number density :$(d)", str_Ar_number_den,   "  1 / m³")
    println(io, "$(b)      He number density :$(d)", str_He_number_den,   "  1 / m³")
    print(io,   "$(b)      H  number density :$(d)", str_H_number_den,    "  1 / m³")

    return nothing
end
