Exponential Atmospheric Model
=============================

```@meta
CurrentModule = SatelliteToolboxAtmosphericModels
DocTestSetup = quote
    using SatelliteToolboxAtmosphericModels
end
```

This model assumes we can compute the atmospheric density by:

```math
\rho(h) = \rho_0 \cdot exp \left\lbrace - \frac{h - h_0}{H} \right\rbrace~,
```

where ``\rho_0``, ``h_0``, and ``H`` are parameters obtained from tables. Reference [1]
provides a discretization of those parameters based on the selected height ``h`` that was
obtained after evaluation of some accurate models.

In this package, we can compute the model using the following function:

```julia
function AtmosphericModels.exponential(h::T) where T<:Number
```

where `h` is the desired height [m].

!!! warning
    Notice that this model does not consider important effects such as the Sun activity, the
    geomagnetic activity, the local time at the desired location, and others. Hence,
    although this can be used for fast evaluations, the accuracy is not good.

## Examples

```jldoctest
julia> AtmosphericModels.exponential(700e3)
3.614e-14
```

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
