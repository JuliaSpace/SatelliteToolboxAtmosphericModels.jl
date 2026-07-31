SatelliteToolboxAtmosphericModels.jl Changelog
==============================================

Version 2.0.0
-------------

- ![BREAKING][badge-breaking] The keyword argument `roots_container` of `jr1971` was
  removed. The model now computes the roots of the quartic polynomial using a closed-form
  algorithm that does not allocate. For the same reason, the dependency PolynomialRoots.jl
  and the ForwardDiff.jl / ImplicitDifferentiation.jl extension were removed since the new
  algorithm is natively compatible with automatic differentiation.
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
