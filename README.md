<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxTransformations.jl"><br>
  <small><i>This package is part of the <a href="https://github.com/JuliaSpace/SatelliteToolbox.jl">SatelliteToolbox.jl</a> ecosystem.</i></small>
</p>

SatelliteToolboxAtmosphericModels.jl
====================================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/branch/main/graph/badge.svg?token=oQOhGnQmdG)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxAtmosphericModels.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package implements atmospheric models for the **SatelliteToolbox.jl** ecosystem.
Currently, the following models are available:

- Exponential atmospheric model;
- Jacchia-Roberts 1971;
- [Jacchia-Bowman 2008](http://sol.spacenvironment.net/jb2008/); and
- [NRLMSISE-00](https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php).

## Installation

```julia
julia> using Pkg
julia> Pkg.add("SatelliteToolboxAtmosphericModels")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolboxAtmosphericModels.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolboxAtmosphericModels.jl/stable
