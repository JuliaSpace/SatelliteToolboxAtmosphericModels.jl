# SatelliteToolboxAtmosphericModels.jl

This package implements atmospheric models for the **SatelliteToolbox.jl** ecosystem.
Currently, the following models are available:

- Exponential atmospheric model according to [1];
- Harris-Priester [2] and modified Harris-Priester [3];
- Jacchia 1977 [4];
- Jacchia-Roberts 1971;
- [Jacchia-Bowman 2008](http://sol.spacenvironment.net/jb2008/); and
- [NRLMSISE-00](https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php).

## Installation

```julia
julia> using Pkg
julia> Pkg.add("SatelliteToolboxAtmosphericModels")
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
- **[2]** **Harris, I., Priester, W** (1962). *Time-dependent structure of the upper
  atmosphere*. **Journal of the Atmospheric Sciences**, 19(4), pp. 286-301.
- **[3]** **Hatten, N., & Russell, R. P. (2017)**. *A smooth and robust Harris-Priester
  atmospheric density model for low Earth orbit applications*. **Advances in Space
  Research**, 59(2), 571-586.
- **[4]** **Jacchia, L. G** (1977). *Thermospheric temperature, density and composition:
  New models*. **SAO Special Report #375**.
