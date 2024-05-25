# SatelliteToolboxAtmosphericModels.jl

This package implements atmospheric models for the **SatelliteToolbox.jl** ecosystem.
Currently, the following models are available:

- Exponential atmospheric model according to [1];
- Jacchia-Roberts 1971;
- [Jacchia-Bowman 2008](http://sol.spacenvironment.net/jb2008/); and
- [NRLMSISE-00](https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php).

## Installation

```julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxAtmosphericModels")
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
