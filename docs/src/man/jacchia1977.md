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
validated against its outputs. Notice that the barometric and diffusion equations are
integrated using the composite Gauss-Legendre quadrature, which is more accurate and much
faster than the Boole rule with a fine step used by the reference code.

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

The following keywords are available:

- `geomagnetic_profile::Val`: Profile of the geomagnetic variation of the temperature,
    only available in the `Val(:sr375)` variant. If it is `Val(:constant)`, the entire
    temperature profile is increased by the geomagnetic variation of the exospheric
    temperature, as in the reference implementation [2]. If it is `Val(:tanh)`, the
    increase is weighted by the altitude-dependent profile of eq. 32 of [1].
    (**Default** = `Val(:constant)`)
- `variant::Val`: Assembly of the dynamic model, `Val(:sr375)` or `Val(:stela)` (see
    below).
    (**Default** = `Val(:sr375)`)

If we omit all space indices, the system tries to obtain them automatically for the
selected day `jd` or `instant` using the prescriptions in [1]. However, the indices must be
already initialized using the function `SpaceIndices.init()`. The Gaussian-weighted mean
uses a window of ±213 days truncated at the available data span.

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
- `H_number_density::T`: Number density of H [1 / m³], which is 0 up to 140 km, where the
  model does not include the hydrogen.

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

## STELA Variant

The CNES tools [STELA](https://www.connectbycnes.fr/en/stela) and
[PATRIUS](https://github.com/CNES/patrius) implement a simplified variant of the Jacchia
1977 model [3]: the static model, precomputed as a lookup table, is evaluated at a single
local exospheric temperature computed with the hydrogen phase angle (-60°) and a fixed
diurnal exponent of 3, the geomagnetic variation of the exospheric temperature is weighted
by the altitude profile of eq. 32 of [1], only the semiannual variation is applied to the
number densities, and the daily and averaged fluxes are swapped in eq. 20 of [1]. This
assembly produces total densities a few percent higher on average than the report
formulation, which directly affects, for example, orbital decay analyses.

The keyword `variant` selects the assembly: `Val(:sr375)` (default) follows the report
[1], whereas `Val(:stela)` follows the CNES tools, allowing the reproduction of analyses
performed with them:

```@repl jacchia1977
AtmosphericModels.jacchia1977(
    DateTime("2018-06-19T18:35:00"),
    deg2rad(-22),
    deg2rad(-45),
    700e3,
    100,
    100,
    3;
    variant = Val(:stela)
)
```

Our implementation of the variant evaluates the static model directly instead of
interpolating the lookup table. It matches the total density of the reference Java
implementation [3] within 0.3 %, which is the bilinear interpolation error of the table
grid, and reproduces the decay time of a 500 km sun-synchronous satellite computed by
STELA within 0.5 %. Notice that, in this variant, the field `exospheric_temperature` of
the output holds the local exospheric temperature used to evaluate the static model
instead of the mean exospheric temperature `T½`, and the keyword `geomagnetic_profile` is
not supported.

## References

- **[1]** **Jacchia, L. G** (1977). *Thermospheric temperature, density and composition:
  New models*. **SAO Special Report #375**.
- **[2]** **de Matos, B. S., Carrara, V** (1985-1987). *Fortran implementation of the
  Jacchia 1977 model*. **INPE**, São José dos Campos, BR. Available in
  [INPE-atmosphere-models](https://github.com/jacobwilliams/INPE-atmosphere-models).
- **[3]** **CNES** (2025). *Java implementation of the Jacchia 1977 model used by the
  tools STELA and PATRIUS* (class fr.cnes.sirius.patrius.stela.forces.atmospheres.
  Jacchia77, PATRIUS 4.16). Available in [CNES/patrius](https://github.com/CNES/patrius).
