# NRLMSISE-00

```@meta
CurrentModule = SatelliteToolboxAtmosphericModels
```

```@setup nrlmsise00
using SatelliteToolboxAtmosphericModels
```

The NRLMSISE-00 empirical atmosphere model was developed by Mike Picone, Alan Hedin, and
Doug Drob based on the MSISE90 model:

> **Picone, J. M., Hedin, A. E., Drob, D. P., Aikin, A. C (2002)**. *NRLMSISE-00 empirical
> model of the atmosphere: Statistical comparisons and scientific issues*. **Journal of
> Geophysical Research: Space Physics**, Vol. 107 (A12), p. SIA 15-1 -- SIA 15-16, DOI:
> 10.1029/2002JA009430.
 
In this package, we can compute this model using the following functions:

```julia
AtmosphericModels.nrlmsise00(instant::DateTime, h::Number, ϕ_gd::Number, λ::Number[, F10ₐ::Number, F10::Number, ap::Union{Number, AbstractVector}]; kwargs...) -> Nrlmsise00Output{Float64}
AtmosphericModels.nrlmsise00(jd::Number, h::Number, ϕ_gd::Number, λ::Number[, F10ₐ::Number, F10::Number, ap::Union{Number, AbstractVector}]; kwargs...) -> Nrlmsise00Output{Float64}
```

where

- `instant::DateTime`: Instant to compute the model represent using `DateTime`.
- `jd::Number`: Julian day to compute the model.
- `h::Number`: Altitude [m].
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `F10ₐ::Number`: 10.7-cm averaged solar flux, 90-day centered on input time [sfu].
- `F10::Number`: 10.7-cm solar flux [sfu].
- `ap::Union{Number, AbstractVector}`: Magnetic index, see the section [AP](@ref) for more
  information.
  
The following keywords are available:

- `flags::Nrlmsise00Flags`: A list of flags to configure the model. For more information,
    see [`AtmosphericModels.Nrlmsise00Flags`](@ref).
    (**Default** = `Nrlmsise00Flags()`)
- `include_anomalous_oxygen::Bool`: If `true`, the anomalous oxygen density will be included
    in the total density computation.
    (**Default** = `true`)
- `P::Union{Nothing, Matrix}`: If the user passes a matrix with dimensions equal to or
    greater than 8 × 4, it will be used when computing the Legendre associated functions,
    reducing allocations and improving the performance. If it is `nothing`, the matrix is
    allocated inside the function.
    (**Default** `nothing`)

If we omit all space indices, the system tries to obtain them automatically for the selected
day `jd` or `instant`. However, the indices must be already initialized using the function
`SpaceIndices.init()`.

These functions return an object of type `Nrlmsise00Output{Float64}` that contains the
following fields:

- `total_density::T`: Total mass density [kg / m³].
- `temperature`: Temperature at the selected altitude [K].
- `exospheric_temperature`: Exospheric temperature [K].
- `N_number_density`: Nitrogen number density [1 / m³].
- `N2_number_density`: N₂ number density [1 / m³].
- `O_number_density`: Oxygen number density [1 / m³].
- `aO_number_density`: Anomalous Oxygen number density [1 / m³].
- `O2_number_density`: O₂ number density [1 / m³].
- `H_number_density`: Hydrogen number density [1 / m³].
- `He_number_density`: Helium number density [1 / m³].
- `Ar_number_density`: Argon number density [1 / m³].

## AP

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

## Examples

```@repl nrlmsise00
AtmosphericModels.nrlmsise00(
    DateTime("2018-06-19T18:35:00"),
    700e3,
    deg2rad(-22),
    deg2rad(-45),
    73.5,
    79,
    5.13
)
```

```@repl nrlmsise00
SpaceIndices.init()

AtmosphericModels.nrlmsise00(DateTime("2018-06-19T18:35:00"), 700e3, deg2rad(-22), deg2rad(-45))
```

If we use the automatic space index fetching mechanism, it is possible to obtain the fetched
values by turning on the debugging logs according to the [Julia
documentation](https://docs.julialang.org/en/v1/stdlib/Logging/):

```@repl nrlmsise00
using Logging

with_logger(ConsoleLogger(stderr, Logging.Debug)) do
    AtmosphericModels.nrlmsise00(DateTime("2018-06-19T18:35:00"), 700e3, deg2rad(-22), deg2rad(-45))
end
```
