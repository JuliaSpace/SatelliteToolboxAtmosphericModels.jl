SatelliteToolboxAtmosphericModels.jl Changelog
==============================================

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
  allocations, which increase the performance by roughly 20%.

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

