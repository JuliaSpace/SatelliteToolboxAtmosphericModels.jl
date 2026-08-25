# Jacchia 1977

```@meta
CurrentModule = SatelliteToolboxAtmosphericModels
```

```@setup jacchia1977
using SatelliteToolboxAtmosphericModels
```

The Jacchia 1977 model [1] is the last of the static models developed by L. G. Jacchia. It
models the atmosphere between 90 km and 2000 km using per-species diffusion equations for
N₂, O₂, O, Ar, He, and H, together with empirical corrections for the diurnal, geomagnetic,
seasonal-latitudinal, and semiannual variations. Its signature feature is the per-species
pseudo exospheric temperature that shifts the diurnal density maximum according to the
molecular mass of each species.

Unlike the Jacchia-Roberts 1971 model, the Jacchia 1977 model does not have a closed-form
solution. Hence, this package numerically integrates the barometric and diffusion
equations, making this model considerably slower than [`AtmosphericModels.jr1971`](@ref).
The implementation follows the reference Fortran code developed at INPE [2] and is
validated against its outputs.

In this package, we can evaluate the model using the following functions:

```julia
AtmosphericModels.jacchia1977(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number[, F10::Number, F10ₐ::Number, Kp::Number]; kwargs...) -> Jacchia1977Output
AtmosphericModels.jacchia1977(jd::Number, ϕ_gd::Number, λ::Number, h::Number[, F10::Number, F10ₐ::Number, Kp::Number]; kwargs...) -> Jacchia1977Output
```

where:

- `jd::Number`: Julian day to compute the model.
- `instant::DateTime`: Instant to compute the model represented using `DateTime`.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `h::Number`: Altitude [m], which must be between 90 km and 2000 km.
- `F10::Number`: 10.7-cm solar flux [sfu], evaluated with the lag prescribed in [1] (from
  0.9 day to 1.6 days, depending on the local solar time).
- `F10ₐ::Number`: 10.7-cm averaged solar flux, Gaussian-weighted mean centered on the input
  time with a standard width of 71 days [sfu].
- `Kp::Number`: Kp geomagnetic index, delayed by 0.1 day to 0.3 day depending on the
  geomagnetic latitude, as prescribed in [1].

If we omit all space indices, the system tries to obtain them automatically for the
selected day `jd` or `instant` using the prescriptions in [1]. However, the indices must be
already initialized using the function `SpaceIndices.init()`, and the Gaussian-weighted
mean requires the indices to be available up to 213 days after the input time.

These functions return an object of type `Jacchia1977Output` that contains the following
fields:

- `total_density::T`: Total atmospheric density [kg / m³].
- `temperature::T`: Local temperature at the selected position [K], including the diurnal
  and geomagnetic variations.
- `exospheric_temperature::T`: Mean exospheric temperature `T½` above the selected position
  [K].
- `N2_number_density::T`: Number density of N₂ [1 / m³].
- `O2_number_density::T`: Number density of O₂ [1 / m³].
- `O_number_density::T`: Number density of O [1 / m³].
- `Ar_number_density::T`: Number density of Ar [1 / m³].
- `He_number_density::T`: Number density of He [1 / m³].
- `H_number_density::T`: Number density of H [1 / m³].

## Examples

```@repl jacchia1977
AtmosphericModels.jacchia1977(
    DateTime("2018-06-19T18:35:00"),
    deg2rad(-22),
    deg2rad(-45),
    700e3,
    100,
    100,
    3
)
```

```@repl jacchia1977
SpaceIndices.init()

AtmosphericModels.jacchia1977(
    DateTime("2023-01-01T10:00:00"),
    deg2rad(-22),
    deg2rad(-45),
    700e3
)
```

If we use the automatic space index fetching mechanism, it is possible to obtain the
fetched values by turning on the debugging logs according to the [Julia
documentation](https://docs.julialang.org/en/v1/stdlib/Logging/):

```@repl jacchia1977
using Logging

with_logger(ConsoleLogger(stderr, Logging.Debug)) do
    AtmosphericModels.jacchia1977(
        DateTime("2023-01-01T10:00:00"),
        deg2rad(-22),
        deg2rad(-45),
        700e3
    )
end
```

## References

- **[1]** **Jacchia, L. G** (1977). *Thermospheric temperature, density and composition:
  New models*. **SAO Special Report #375**.
- **[2]** **de Matos, B. S., Carrara, V** (1985-1987). *Fortran implementation of the
  Jacchia 1977 model*. **INPE**, São José dos Campos, BR. Available in
  [INPE-atmosphere-models](https://github.com/jacobwilliams/INPE-atmosphere-models).
