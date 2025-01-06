SatelliteToolboxAtmosphericModels.jl Changelog
==============================================

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

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

[gh-pr-3]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/pull/3
