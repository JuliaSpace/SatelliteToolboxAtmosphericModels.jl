SatelliteToolboxAtmosphericModels.jl Changelog
==============================================

Version 2.0.1
-------------

- ![Bugfix][badge-bugfix] Bump SatelliteToolboxBase.jl to v2.


Version 2.0.0
-------------

- ![BREAKING][badge-breaking] The keyword `n` (cosine exponent of the diurnal bulge) of both
  Harris-Priester models now accepts any `Number` in the interval `[2, 7]`, throwing an
  `ArgumentError` otherwise. Previously, the classic model required an `Int` in `[2, 6]`,
  and the modified model accepted any value without validation. The exponent no longer
  participates in the output type promotion of the modified model.
- ![BREAKING][badge-breaking] The field `turbo_scale_height` of `Nrlmsise00Flags` was
  removed since it had no effect on the model (it is also unused in the reference
  implementation).
- ![BREAKING][badge-breaking] All the models now validate the altitude at their entry
  points, throwing an `ArgumentError` outside their validity ranges: `[0, ∞)` for the
  exponential and NRLMSISE-00 models, `[100, 1000]` km for the Harris-Priester models (the
  range of the density profile for the classic model), `[90, 3000]` km for JR1971 and
  JB2008, and `[90, 2000]` km for Jacchia 1977. Previously, the classic Harris-Priester
  model returned zero above 1000 km, the modified Harris-Priester model silently
  extrapolated the density outside its range (returning, e.g., 7.5 kg / m³ at the sea
  level), and JR1971 and JB2008 had no upper bound. The classic Harris-Priester model also
  validates the density profile passed with the keyword `alt_ρ`.
- ![BREAKING][badge-breaking] The deprecated keyword `roots_container` of `jr1971`, which
  had no effect since v1.3.0, was removed.
- ![BREAKING][badge-breaking] The keyword `verbose` was removed from all the models. The
  debug messages related to the automatic space index fetching are now always emitted
  through the logging system, which can be configured to show or hide them.
- ![Enhancement][badge-enhancement] The output structures of the models are now subtypes of
  the new abstract type `AbstractAtmosphericModelOutput`, and their `show` methods share a
  single implementation. The dependency on Crayons.jl was removed.
- ![Enhancement][badge-enhancement] The NRLMSISE-00 model computes the cubic spline of the
  lower thermosphere temperature profile once per evaluation instead of once per species,
  making the model about 13 % faster below 72.5 km. The Legendre functions are now stored in
  a static matrix and the internal helpers were simplified. The results are unchanged.
- ![Enhancement][badge-enhancement] The JB2008 model evaluates the polynomials of the local
  solar time and latitude correction of the exospheric temperature once instead of up to
  three times, validates the exospheric temperature once instead of at every quadrature
  node, and uses cheaper expressions for the fourth root of the flux ratio and for the
  number of integration steps. The model is about 3 % faster, and the results are unchanged
  to the round-off.
- ![Enhancement][badge-enhancement] The JR1971 model is about 35 % faster above 125 km and
  10 % faster below because the powers of the temperature ratios are evaluated as
  exponentials of shared logarithms, base-10 powers use `exp10`, and the output structure is
  built by a single helper. The results are unchanged to the round-off.
- ![Enhancement][badge-enhancement] The diurnal variation of the Jacchia 1977 model no
  longer integrates the hydrogen for the five static model evaluations that only use the
  heavy species, and the quiet temperature profile is selected by dispatch instead of a
  test at every quadrature node. Together with the new quadrature, the model is about 12
  times faster at 110 km, 5 times faster at 300 km, and 2.4 times faster at 1500 km than
  the previous version.
- ![Enhancement][badge-enhancement] The static model of Jacchia 1977 now integrates the
  barometric and diffusion equations using the composite 8-point Gauss-Legendre quadrature
  with the gravity evaluated at the nodes instead of the Boole rule with a fine step and
  the gravity frozen at the center of each panel, as in the reference Fortran
  implementation. The new scheme reproduces a converged integration to about 1e-7 in the
  base-10 logarithm of the number densities of the heavy species, whereas the previous one
  deviated by up to 2e-4. The hydrogen flux term is now integrated with the Heun method
  (trapezoidal rule for the ratio between the heavy species and the hydrogen number
  densities), reproducing a converged integration to about 7e-3 in the base-10 logarithm
  of the hydrogen number density with 10 km panels. The model is 2.6 times faster at
  300 km, and the total density changes by less than 0.01 % with respect to the previous
  version.
- ![Enhancement][badge-enhancement] The modified Harris-Priester model evaluates the diurnal
  factor first and skips the maximum density profile on the night side, and the minimum and
  maximum profiles share a single implementation. The model is about 24 % faster, and the
  results are unchanged.
- ![Enhancement][badge-enhancement] The tables of the exponential and Harris-Priester models
  are now static arrays instead of mutable global vectors and matrices, and the density
  profile of the classic Harris-Priester model is searched with a row-indexed binary search
  that works for any matrix type without allocating.
- ![Bugfix][badge-bugfix] The keyword `P` of `nrlmsise00` silently copied the pre-allocated
  matrix (allocating on every call) when its type was not `Matrix{T}`, where `T` is the
  promoted input type, and truncated the Legendre functions to the element type of `P`
  when it was narrower. The keyword now accepts any `AbstractMatrix` without copying and
  throws an `ArgumentError` if its element type differs from the promoted input type. The
  dependency on LinearAlgebra.jl was removed.
- ![Bugfix][badge-bugfix] The NRLMSISE-00 model threw an `InexactError` when all the inputs
  were integers. The output element type is now the promotion of the input types converted
  to a floating-point type, as in the other models.
- ![Bugfix][badge-bugfix] The NRLMSISE-00 model converted the input latitude and longitude
  to degrees using the truncated internal constant of the reference implementation
  (`1.74533e-2`) instead of the exact conversion, introducing a relative error of `4e-7` in
  the angles. The inputs are now converted exactly, changing the results by a negligible
  amount (about `1e-7` relative).
- ![Bugfix][badge-bugfix] The NRLMSISE-00 model now throws an `ArgumentError` if the vector
  `ap` does not have 7 elements, instead of a `BoundsError` from a private function.
- ![Bugfix][badge-bugfix] The NRLMSISE-00 model returned `NaN` for all densities below the
  mesopause (72.5 km) when the flag `departures_from_eq` was `false` because the mixed N₂
  density used to blend the thermospheric and lower atmosphere profiles was not computed.
- ![Bugfix][badge-bugfix] The number densities of the species returned by the JR1971 model
  between 90 km and 100 km were computed from the total density using the constituent
  fractions of the model without the molecular mass factor used above 100 km. Hence, the
  species mass densities summed to 104.8 % of the total density, and the number density of,
  e.g., atomic oxygen was 81 % higher than the consistent value. The species are now
  computed as in the region between 100 km and 125 km. The total density is unchanged.
- ![Bugfix][badge-bugfix] The Jacchia 1977 model returned a hydrogen number density of about
  1 / m³ up to 140 km, where the model does not include the hydrogen, due to an internal
  placeholder. It now returns 0, as the JR1971 model does below 500 km. Additionally, the
  local temperature computed exactly at 90 km now keeps the type of the altitude, fixing
  the derivative with respect to the altitude at that point in automatic differentiation.
- ![Bugfix][badge-bugfix] The JR1971 model returned wrong densities between 90 km and
  100 km. The closed-form solution of the barometric equation presented in the reference
  (and also implemented in GMAT) led to an almost constant density in this region and to a
  discontinuity of a factor of about 6 at 100 km (e.g. 3.07e-6 kg / m³ at 99.999 km and
  5.38e-7 kg / m³ at 100.001 km). The model now integrates the barometric equation
  numerically using the 8-point Gauss-Legendre quadrature, which reproduces a fine
  numerical integration to better than 1e-10 and makes the density continuous at 90 km and
  100 km. The results above 100 km are unchanged.
- ![Bugfix][badge-bugfix] The exponential model returned a `Float64` for `Float32` inputs.
  It now returns the floating-point type of the input.
- ![Bugfix][badge-bugfix] The Zygote.jl rule of `nrlmsise00` now supports the 7-element
  magnetic index vector `ap`. Previously, only the daily index was supported, and Zygote.jl
  failed for the vector input.
- ![Info][badge-info] The allocation tests now verify at runtime that every model is
  allocation-free, including the methods that fetch the space indices automatically, and
  the static checks with AllocCheck.jl cover the methods with explicit space indices.

Version 1.5.0
-------------

- ![Feature][badge-feature] The Jacchia 1977 model now provides the keyword `variant` to
  select the assembly of the dynamic model. The default, `Val(:sr375)`, keeps the
  formulation of the report. The new option, `Val(:stela)`, follows the simplified
  assembly used by the CNES tools STELA and PATRIUS: the static model is evaluated at a
  single local exospheric temperature computed with the hydrogen phase angle and a fixed
  diurnal exponent, the geomagnetic variation of the exospheric temperature is weighted by
  the altitude profile of eq. 32 of the report, only the semiannual variation is applied
  to the number densities, and the daily and averaged fluxes are swapped in eq. 20, as in
  the reference Java implementation. This variant produces total densities a few percent
  higher on average, and it allows reproducing analyses performed with those tools: the
  decay time of a 500 km sun-synchronous satellite computed by STELA is reproduced within
  0.5 %, whereas the report formulation yields a decay time about 7 % longer. The variant
  is validated against the density table and the algorithm of the class
  fr.cnes.sirius.patrius.stela.forces.atmospheres.Jacchia77 of PATRIUS 4.16.

Version 1.4.0
-------------

- ![Bugfix][badge-bugfix] The automatic space index fetching of the JR1971 and Jacchia
  1977 models used the observed F10.7 flux, whereas the Jacchia models were fitted with
  the flux adjusted to 1 AU, as documented in the reference Fortran implementation of the
  Jacchia 1977 model. The fetching now uses the adjusted flux. The two series differ by up
  to ±3.4 %, which maps to several percent in the auto-fetched density near the perihelion
  and aphelion. Users passing the indices manually are not affected.
- ![Bugfix][badge-bugfix] The automatic space index fetching of the NRLMSISE-00 model used
  the F10.7 flux adjusted to 1 AU, whereas the model documentation explicitly requires the
  observed flux at the actual distance of the Earth from the Sun (as also stated in the
  docstring). The fetching now uses the observed flux. The two series differ by up to
  ±3.4 %, which maps to a spurious annual signature of up to ±7 % in the auto-fetched
  density near the perihelion and aphelion. Users passing the indices manually are not
  affected.
- ![Enhancement][badge-enhancement] Several type instabilities related to `Float64`
  literals were removed from the JR1971, JB2008, and modified Harris-Priester models. In
  particular, the modified Harris-Priester model now returns the promotion of the input
  types (e.g. `Float32` for all-`Float32` inputs) as documented, instead of always
  returning `Float64`. The results for `Float64` inputs are unchanged.
- ![Enhancement][badge-enhancement] The NRLMSISE-00 coefficient tables are now tuples
  instead of heap-allocated vectors. Since all the accesses use literal indices, the bounds
  checks are elided at compile time, making the model about 3 % faster. The results are
  unchanged.
- ![Bugfix][badge-bugfix] The standard library `Dates` is now a declared dependency. It was
  previously reachable only through a re-export of SatelliteToolboxBase.jl, which could
  break silently if that upstream re-export changed.
- ![Bugfix][badge-bugfix] The JR1971, JB2008, and modified Harris-Priester models threw an
  `InexactError` when all the inputs were integers. The output element type is now the
  promotion of the input types converted to a floating-point type. Additionally, the
  Harris-Priester cosine exponent is now validated with an `ArgumentError` instead of an
  `@assert`, which can be disabled by compiler options.
- ![Bugfix][badge-bugfix] The JB2008 model now validates the altitude at the entry point,
  throwing an `ArgumentError` for altitudes below 90 km, and its docstring documents the
  bound. Previously, such altitudes emitted a spurious warning and an error message
  mentioning an altitude of 0 km.
- ![Bugfix][badge-bugfix] The JR1971 model returned `NaN` for all densities when the Sun
  declination was exactly zero (equinox) because the helium seasonal correction contained a
  0 / 0 term. The correction is now written using `sign`, which is finite and returns the
  same values elsewhere.
- ![Bugfix][badge-bugfix] The `temperature` field returned by the Jacchia 1977 model was
  computed from the temperature profile related to the mean exospheric temperature `T½`,
  ignoring the diurnal and geomagnetic variations, as in the reference Fortran
  implementation. It is now the local temperature, computed from the profile related to the
  local quiet exospheric temperature (phase angle of -60°, prescribed for the actual
  temperature in eq. 26 of the report) increased by the geomagnetic variation. For the
  worked example of the report, the returned temperature at high altitudes now matches the
  report values (939.3 K quiet and 1061 K disturbed) to better than 0.5 K. The densities
  are unaffected.
- ![Feature][badge-feature] The Jacchia 1977 model now supports the keyword
  `geomagnetic_profile` to select how the geomagnetic variation of the exospheric
  temperature is applied to the temperature profile. The default, `Val(:constant)`,
  increases the entire profile, as in the reference Fortran implementation and in the
  numerical example of the report. `Val(:tanh)` weights the increase by the
  altitude-dependent profile of eq. 32 of the report, which the report states is required
  at lower heights.
- ![Bugfix][badge-bugfix] The Jacchia 1977 model clamped negative base-10 logarithms of the
  number densities to 0, as in the reference Fortran implementation. Hence, the number
  density of heavily depleted species (e.g. Ar, O₂, and N₂ at high altitudes) was reported
  as 1 / m³ instead of its true, much smaller value. The total density is essentially
  unaffected.
- ![Bugfix][badge-bugfix] The hydrogen integration of the Jacchia 1977 model reproduced two
  inaccuracies of the reference Fortran implementation: below 500 km, the hydrogen was
  integrated over a window displaced by one Boole panel (about 20 km) from the window used
  by the other species; and above 500 km, the flux term of eq. 16 of the report omitted the
  total number density factor. Both were fixed, changing the hydrogen number density by up
  to 6 % (at 150 km) with respect to the previous version. The total density changes by
  less than 0.25 %, and only above approximately 1500 km, where the hydrogen dominates the
  mass.
- ![Bugfix][badge-bugfix] The automatic space index fetching of the Jacchia 1977 model
  selected the daily F10.7 flux of the day before the lagged instant prescribed in eq. 23
  of the report because SpaceIndices.jl shifts the F10.7 lookups by -8 hours to center the
  intervals on the measurement time. The fetching now compensates this shift, returning the
  flux tabulated for the calendar day of the lagged instant.
- ![Bugfix][badge-bugfix] The automatic space index fetching of the Jacchia 1977 model
  required the F10.7 flux to be available up to 213 days after the input time, making the
  model unusable for recent epochs (the fetching threw an `ArgumentError`). The Gaussian
  window is now truncated at the available data span and the weights are renormalized
  accordingly, following the definition of the weighted mean in eq. 21 of the report.
- ![Bugfix][badge-bugfix] The numerical integrator of the Jacchia 1977 model could loop
  forever for reduced-precision inputs (e.g. `Float32`) because the loop termination
  compared the accumulated altitude against a fixed tolerance of 1e-4 km, which is smaller
  than the accumulated rounding error of such types. The loops now iterate over an integer
  panel count. The results for `Float64` inputs are unchanged.

Version 1.3.0
-------------

- ![Deprecation][badge-deprecation] The keyword argument `roots_container` of `jr1971` is
  not used anymore and is kept only for backward compatibility. The model now computes the
  roots of the quartic polynomial using a closed-form algorithm that does not allocate.
  For the same reason, the dependency PolynomialRoots.jl and the ForwardDiff.jl /
  ImplicitDifferentiation.jl extension were removed since the new algorithm is natively
  compatible with automatic differentiation.
- ![Feature][badge-feature] The package now supports the Jacchia 1977 model
  (`AtmosphericModels.jacchia1977`). The implementation numerically integrates the
  barometric and diffusion equations as described in SAO Special Report #375 and is
  validated against the reference Fortran implementation developed at INPE. The automatic
  space index fetching follows the prescriptions in the report: lagged daily F10.7,
  Gaussian-weighted averaged F10.7, and Kp delayed by a geomagnetic latitude dependent
  interval.
- ![Bugfix][badge-bugfix] The nighttime minimum global exospheric temperature in the JR1971
  model used the daily F10.7 flux in the term that requires the 81-day average, leading to
  wrong densities whenever both indices differ.
- ![Bugfix][badge-bugfix] The automatic space index fetching in the JR1971 model computed
  the 3-hour delayed Kp using the slot index of the non-delayed instant. Hence, the delay
  was effectively not applied and, near day boundaries, the selected Kp could be off by
  almost one day.
- ![Bugfix][badge-bugfix] The hydrogen number density of the JB2008 model was approximately
  1 / m³ below 105 km due to an assignment to an unused variable. It now equals the helium
  number density divided by exp(25), as in the reference implementation. The total density
  is essentially unaffected.
- ![Bugfix][badge-bugfix] The Earth radius used in the JB2008 gravity computation had a
  digit typo (6356.776 km instead of 6356.766 km), causing a relative bias of about 1e-6 in
  the barometric integrals.
- ![Bugfix][badge-bugfix] The NRLMSISE-00 helium mixing correction was gated by the N₂
  altitude limit (160 km) instead of the helium limit (200 km), changing the helium and
  total densities between 160 km and 200 km when compared with the reference C
  implementation. Additionally, the N₂ branch now uses the same inclusive comparison as the
  reference.
- ![Bugfix][badge-bugfix] The Harris-Priester model could throw a `DomainError` due to
  rounding errors when the angle between the diurnal bulge apex and the satellite was close
  to 180°.
- ![Bugfix][badge-bugfix] The NRLMSISE-00 model raised a `MethodError` for non-Float64
  inputs (e.g. `Float32`) below the mesopause due to heterogeneous temperature tuples.
- ![Bugfix][badge-bugfix] The Zygote.jl rule for `nrlmsise00` returned a wrong number of
  tangents, failed when a pre-allocated Legendre matrix was provided, and could not handle
  the AP vector input. The rule was fixed, restricted to scalar AP, and now computes the
  full Jacobian in a single ForwardDiff.jl sweep instead of one sweep per output field.
  Additionally, the `_get_doy` rule, previously duplicated in the Zygote.jl and Mooncake.jl
  extensions, now lives in a new ChainRulesCore.jl extension shared by both backends.
- ![Enhancement][badge-enhancement] The automatic space index fetching in the NRLMSISE-00
  model now uses the 81-day centered average of the F10.7 flux stated in the model
  documentation instead of a 90-day non-centered window.
- ![Enhancement][badge-enhancement] The NRLMSISE-00 model no longer mutates the shared
  coefficient vectors when clamping two parameters, making the function thread-safe.
- ![Enhancement][badge-enhancement] `jr1971` is now allocation-free and all models are
  faster due to several improvements: closed-form quartic solver, hoisting of repeated
  subexpressions, integer powers instead of float-exponent powers, binary search in the
  exponential model table lookup, and removal of dynamic dispatches in the modified
  Harris-Priester model.
- ![Enhancement][badge-enhancement] Several type instabilities were removed, improving the
  support for custom number types (e.g. dual numbers used by automatic differentiation).
- ![Enhancement][badge-enhancement] The automatic space index fetching in the modified
  Harris-Priester model now logs the fetched indices using debug messages, as the other
  models.
- ![Info][badge-info] Many typos, grammar errors, and documentation errors were fixed,
  including wrong units in the output structures, wrong return types, and undocumented
  keyword arguments.

Version 1.2.1
-------------

- ![Info][badge-info] We updated the license information and links for the JB2008 model.

Version 1.2.0
-------------

- ![Feature][badge-feature] The package now supports the Harris Priester model. (PR
  [#9][gh-pr-9] and [#10][gh-pr-10])

Version 1.1.1
-------------

- ![Bugfix][badge-bugfix] 𝜏 angle must be in the interval [-π, π] in JR1971 model. (PR
  [#8][gh-pr-8])

Version 1.1.0
-------------

- ![Feature][badge-feature] The package now supports automatic differentiation using
  different backends. (PR [#6][gh-pr-6])
- ![Enhancement][badge-enhancement] Some allocations were removed. (PR [#6][gh-pr-6])

Version 1.0.0
-------------

- ![Enhancement][badge-enhancement] The functions now support automatic differentiation.
  (PR [#3][gh-pr-3])
- ![Info][badge-info] We dropped support for Julia 1.6. This version only supports the
  current Julia version and v1.10 (LTS).
- ![Info][badge-info] This version does not have breaking changes. We bump the version to
  1.0.0 because we now consider the API stable.

Version 0.1.3
-------------

- ![Enhancement][badge-enhancement] Minor source-code updates.
- ![Enhancement][badge-enhancement] Documentation updates.

Version 0.1.2
-------------

- ![Bugfix][badge-bugfix] In certain conditions, the JR1971 model would generate an array
  index that was not an integer. This bug is now fixed.
- ![Feature][badge-feature] NRLMSISE-00 can receive a matrix to call the in-place function
  that computes the Legendre associated functions, reducing the allocations.
- ![Enhancement][badge-enhancement] The internal structure of NRLMSISE-00 model was changed
  to have a type parameter to indicate whether the AP is a vector or number. This approach
  slightly increased the compile time, but reduced one allocation.
- ![Enhancement][badge-enhancement] The user can now call NRLMSISE-00 model without any
  allocations, which increases the performance by roughly 20%.

Version 0.1.1
-------------

- ![Enhancement][badge-enhancement] We updated the dependency compatibility bounds.

Version 0.1.0
-------------

- Initial version.
  - This version was based on the functions in **SatelliteToolbox.jl**.

[badge-breaking]: https://img.shields.io/badge/Breaking-DC2626?style=flat-square
[badge-deprecation]: https://img.shields.io/badge/Deprecation-D97706?style=flat-square
[badge-feature]: https://img.shields.io/badge/Feature-16A34A?style=flat-square
[badge-enhancement]: https://img.shields.io/badge/Enhancement-0284C7?style=flat-square
[badge-bugfix]: https://img.shields.io/badge/Bugfix-DB2777?style=flat-square
[badge-info]: https://img.shields.io/badge/Info-475569?style=flat-square

[gh-pr-3]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/pull/3
[gh-pr-6]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/pull/6
[gh-pr-8]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/pull/8
[gh-pr-9]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/pull/9
[gh-pr-10]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/pull/10
