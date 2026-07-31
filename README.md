<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxAtmosphericModels.jl"><br>
  <small><i>This package is part of the <a href="https://github.com/JuliaSpace/SatelliteToolbox.jl">SatelliteToolbox.jl</a> ecosystem.</i></small>
</p>

# SatelliteToolboxAtmosphericModels.jl

[![CI](https://img.shields.io/github/actions/workflow/status/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/ci.yml?style=flat-square&logo=githubactions&logoColor=white&labelColor=475569&label=CI)](https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/github/JuliaSpace/SatelliteToolboxAtmosphericModels.jl?token=oQOhGnQmdG&style=flat-square&logo=codecov&logoColor=white&labelColor=475569)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxAtmosphericModels.jl)
[![docs-stable](https://img.shields.io/badge/docs-stable-16A34A?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-stable-url]
[![docs-dev](https://img.shields.io/badge/docs-dev-D97706?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495D1?style=flat-square&logo=julia&logoColor=white&labelColor=475569)](https://github.com/invenia/BlueStyle)
[![License](https://img.shields.io/github/license/JuliaSpace/SatelliteToolboxAtmosphericModels.jl?style=flat-square&logo=readme&logoColor=white&labelColor=475569&color=0284C7)](https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/blob/main/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.10644917-DB2777?style=flat-square&logo=doi&logoColor=white&labelColor=475569)](https://zenodo.org/doi/10.5281/zenodo.10644917)

This package implements atmospheric models for the **SatelliteToolbox.jl** ecosystem.
Currently, the following models are available:

- Exponential atmospheric model;
- Harris-Priester;
- Modified Harris-Priester;
- Jacchia 1977;
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
