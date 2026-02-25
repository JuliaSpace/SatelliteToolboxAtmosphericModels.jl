# Harris-Priester

```@meta
CurrentModule = SatelliteToolboxAtmosphericModels
```

```@setup harrispriester
using SatelliteToolboxAtmosphericModels
```

The Harris-Priester model is a simple, computationally efficient atmospheric density model
for altitudes between 100 km and 1000 km. It models the diurnal density variation using a
cosine power law and exponential interpolation between tabulated minimum and maximum density
profiles.

This package provides two variants of the model: the classic Harris-Priester and a modified
version published in:

> **Hatten, N., & Russell, R. P. (2017)**. *A smooth and robust Harris-Priester atmospheric
> density model for low Earth orbit applications*. **Advances in Space Research**, 59(2),
> 571-586.

The modified variant uses cubic polynomial fits on the 81-day centered F10.7 solar flux
average to model density variations across solar activity levels, ensuring continuous first
derivatives and eliminating singularities.

## Harris-Priester (Classic)

In this package, we can evaluate the classic model using the following functions:

```julia
AtmosphericModels.harrispriester(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number; kwargs...) -> Float64
AtmosphericModels.harrispriester(jd::Number, ϕ_gd::Number, λ::Number, h::Number; kwargs...) -> Float64
```

where:

- `jd::Number`: Julian day to compute the model.
- `instant::DateTime`: Instant to compute the model represent using `DateTime`.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `h::Number`: Altitude [m].

The following keywords are available:

- `n::Int`: Cosine exponent for the diurnal bulge (`2` ≤ `n` ≤ `6`). A value of `2`
    produces a smooth day-night transition; `6` produces a sharp transition.
    (**Default** = `4`)
- `alt_ρ::AbstractMatrix`: Custom density profile table with columns `[altitude [m],
    ρ_min [kg/m³], ρ_max [kg/m³]]`.
    (**Default** = built-in table for mean solar activity)

These functions return the atmospheric density [kg/m³].

### Examples

```@repl harrispriester
AtmosphericModels.harrispriester(
    DateTime("2018-06-19T18:35:00"),
    deg2rad(-22),
    deg2rad(-45),
    700e3
)
```

```@repl harrispriester
AtmosphericModels.harrispriester(
    DateTime("2018-06-19T18:35:00"),
    deg2rad(-22),
    deg2rad(-45),
    700e3;
    n = 2
)
```

## Modified Harris-Priester

In this package, we can evaluate the modified model using the following functions:

```julia
AtmosphericModels.harrispriester_modified(instant::DateTime, ϕ_gd::Number, λ::Number, h::Number[, F10ₐ::Number]; n::Number = 4) -> Float64
AtmosphericModels.harrispriester_modified(jd::Number, ϕ_gd::Number, λ::Number, h::Number[, F10ₐ::Number]; n::Number = 4) -> Float64
```

where:

- `jd::Number`: Julian day to compute the model.
- `instant::DateTime`: Instant to compute the model represent using `DateTime`.
- `ϕ_gd::Number`: Geodetic latitude [rad].
- `λ::Number`: Longitude [rad].
- `h::Number`: Altitude [m].
- `F10ₐ::Number`: 81-day centered average of the F10.7 solar flux index [sfu].
- `n::Number`: Cosine exponent for the diurnal bulge. The original Fortran implementation
    computes this from orbital inclination as `n = 2.001 + 4sin²(i)`, yielding `n ≈ 2` for
    equatorial orbits and `n ≈ 6` for polar orbits.
    (**Default** = `4`)

If we omit all space indices, the system tries to obtain them automatically for the selected
day `jd` or `instant`. However, the indices must be already initialized using the function
`SpaceIndices.init()`.

These functions return the atmospheric density [kg/m³].

### Examples

```@repl harrispriester
AtmosphericModels.harrispriester_modified(
    DateTime("2018-06-19T18:35:00"),
    deg2rad(-22),
    deg2rad(-45),
    700e3,
    79.0
)
```

```@repl harrispriester
SpaceIndices.init()

AtmosphericModels.harrispriester_modified(
    DateTime("2022-06-19T18:35:00"),
    deg2rad(-22),
    deg2rad(-45),
    700e3
)
```

If we use the automatic space index fetching mechanism, it is possible to obtain the fetched
values by turning on the debugging logs according to the [Julia
documentation](https://docs.julialang.org/en/v1/stdlib/Logging/):

```@repl harrispriester
using Logging

with_logger(ConsoleLogger(stderr, Logging.Debug)) do
    AtmosphericModels.harrispriester_modified(
        DateTime("2022-06-19T18:35:00"),
        deg2rad(-22),
        deg2rad(-45),
        700e3
    )
end
```
