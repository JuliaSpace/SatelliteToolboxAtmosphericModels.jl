## Description #############################################################################
#
# Shared functions to show the results of the atmospheric models.
#
############################################################################################

# ANSI escape sequences to print bold text and to reset the text style.
const _B = "\e[1m"
const _D = "\e[0m"

# Formats used to print the values in the `text/plain` representation of the outputs.
const _SHOW_FORMAT_DENSITY     = Printf.Format("%15g")
const _SHOW_FORMAT_TEMPERATURE = Printf.Format("%15.2f")

"""
    _model_name(::Type{<:AbstractAtmosphericModelOutput}) -> String

Return the short name of the atmospheric model related to the output type (e.g.
`"JR1971"`), used by the compact `show` method.
"""
function _model_name end

"""
    _model_description(::Type{<:AbstractAtmosphericModelOutput}) -> String

Return the full name of the atmospheric model related to the output type (e.g.
`"Jacchia-Roberts 1971"`), used by the `text/plain` `show` method.
"""
function _model_description end

"""
    _show_fields(::Type{<:AbstractAtmosphericModelOutput}) -> Tuple

Return a tuple describing the fields shown by the `text/plain` `show` method of the output
type. Each element is a tuple `(label, field, format, unit)`, where `label::String` is the
text printed before the value, `field::Symbol` is the field name, `format::Printf.Format`
is the format used to print the value, and `unit::String` is the unit printed after the
value.
"""
function _show_fields end

"""
    show(io::IO, out::AbstractAtmosphericModelOutput) -> Nothing

Print the compact representation of `out`, containing the model name and the total density
[kg / m³].
"""
function show(io::IO, out::T) where {T <: AbstractAtmosphericModelOutput}
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? _B : ""
    d = color ? _D : ""

    print(
        io,
        b,
        _model_name(T),
        " output",
        d,
        " (ρ = ",
        @sprintf("%g", out.total_density),
        " kg / m³)",
    )

    return nothing
end

"""
    show(io::IO, ::MIME"text/plain", out::AbstractAtmosphericModelOutput) -> Nothing

Print the multi-line representation of `out`, containing one line per field returned by
[`_show_fields`](@ref) with the values right-aligned and followed by their units.
"""
function show(
    io::IO, ::MIME"text/plain", out::T
) where {T <: AbstractAtmosphericModelOutput}
    # Check for color support in the `io`.
    color = get(io, :color, false)
    b = color ? _B : ""
    d = color ? _D : ""

    fields = _show_fields(T)

    # The labels are right-aligned to the longest one plus two leading spaces.
    label_width = maximum(textwidth(f[1]) for f in fields) + 2

    println(io, _model_description(T), " Atmospheric Model Result:")

    for (k, (label, field, format, unit)) in enumerate(fields)
        value = Printf.format(format, getfield(out, field))

        print(io, b, lpad(label, label_width), " :", d, value, "  ", unit)

        # The last line must not end with a newline.
        (k != length(fields)) && println(io)
    end

    return nothing
end
