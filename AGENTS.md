# Repository Guide

SatelliteToolboxAtmosphericModels.jl implements atmospheric density models for the SatelliteToolbox.jl ecosystem: Exponential, Harris-Priester (classic and modified), Jacchia 1977, Jacchia-Roberts 1971, Jacchia-Bowman 2008, and NRLMSISE-00.

## Package Structure

- Single package. Requires Julia 1.10 or newer (`[compat] julia = "1.10"`).
- `src/SatelliteToolboxAtmosphericModels.jl` is the entrypoint. It re-exports SpaceIndices.jl, exports the submodule `AtmosphericModels`, and includes `src/AtmosphericModels.jl`, which defines that submodule and holds every `using` and `include`.
- Include order inside `AtmosphericModels`: `types.jl`, `utils.jl`, then one subdirectory per model (`exponential/`, `harrispriester/`, `jacchia1977/`, `jr1971/`, `jb2008/`, `nrlmsise00/`, each with `constants.jl`, an optional `types.jl`, and the model file), and `show.jl` last. New code can only reference symbols included before it.
- Every model is accessed as `AtmosphericModels.<model>(...)`. `export` statements sit at the top of the file that defines the symbol, never in a centralized list.
- Package extensions in `ext/`, declared in `[weakdeps]`/`[extensions]`: `SatelliteToolboxAtmosphericModelsChainRulesCoreExt` (trigger: ChainRulesCore) defines the shared reverse-mode rules; `...MooncakeExt` (triggers: Mooncake + ChainRulesCore) and `...ZygoteExt` (triggers: Zygote + ForwardDiff) build on it. The Zygote extension differentiates NRLMSISE-00 through a ForwardDiff Jacobian because the model mutates. ForwardDiff works without any extension.
- Tests are wired in `test/runtests.jl`. The six model files (`test/exponential.jl`, `test/harrispriester.jl`, `test/jacchia1977/jacchia1977.jl`, `test/jr1971.jl`, `test/jb2008/jb2008.jl`, `test/nrlmsise00.jl`) always run. `test/performance.jl` (Aqua, JET, AllocCheck) and `test/extensions.jl` (ForwardDiff, Zygote, Mooncake) run only on non-prerelease Julia. JET and AllocCheck checks are additionally skipped on Julia 1.12 or newer.
- `test/jacchia1977/*.dat` and `test/jb2008/*.DAT` are reference outputs of the original Fortran implementations. Treat them as immutable fixtures.
- Most test files call `SpaceIndices.init()`, which downloads the space index files on first use, so the suite needs internet access the first time it runs on a machine.
- Test-only dependencies are declared in `[extras]` + `[targets]` of `Project.toml`. There is no `test/Project.toml`.
- `Manifest.toml` is gitignored and never committed. It is tied to the Julia version that generated it: a Manifest left behind by another Julia version can make stdlib precompilation fail (observed with Julia 1.13 on a Manifest from 1.10). If `using Test` or `using Pkg` fails to precompile under `--project=.`, delete `Manifest.toml` and instantiate again.
- There is no `deps/build.jl`, so `Pkg.build()` does nothing here.

## Commands

- Instantiate: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Full test suite (reproduces CI): `julia --project=. -e 'using Pkg; Pkg.test()'`. The first run precompiles for minutes while printing little; use generous timeouts.
- Focused model test file: `julia --project=. -e 'using SatelliteToolboxAtmosphericModels, SatelliteToolboxBase, Test; include("test/jr1971.jl")'` (same pattern for the other five model files).
- `test/performance.jl` and `test/extensions.jl` need test-only packages that are not available in a plain `--project=.` session; run them through `Pkg.test()`.
- There is no test-name selector: `test/runtests.jl` does not read `ARGS`.
- Build docs: `julia --project=docs docs/make.jl`. First run: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`. Pass `local` as an argument to disable pretty URLs for browsing `docs/build/` (gitignored) offline.
- CI (`.github/workflows/ci.yml`) builds the package (`julia-buildpkg`) and then runs the tests on Julia 1.10 and the latest stable 1.x, on Ubuntu x64, macOS arm64, and Windows x64, uploading coverage to Codecov. `ci-nightly.yml` runs the same matrix on Julia nightly, where the performance and extension tests are skipped by `runtests.jl`. `docs.yml` builds and deploys the documentation on the latest stable Julia. CompatHelper and TagBot are automated.

## Code Style

- Code follows Blue Style with the overrides in `.JuliaFormatter.toml` at the repo root (`style = "blue"`, alignment flags on, `whitespace_in_kwargs = true`). That file is the source of truth; do not edit it.
- Line width is 92 columns with greedy wrapping, 4-space indent, no tabs. Comment headers, section separators, and end-of-line comments are filled to exactly column 92.
- Format: `julia -e 'using JuliaFormatter; format(".")'`. JuliaFormatter is not a project dependency; install it in the default environment first with `julia -e 'using Pkg; Pkg.add("JuliaFormatter")'`. Do not use `--project=.` for it.
- CI does not run a format check. Format manually before committing.
- Every source file starts with a `## Description ####...` header block (plus a `## References ####...` block citing the papers or reference code when applicable), and groups private functions under a `#  Private Functions  #` banner (`src/jr1971/jr1971.jl` shows the full layout; `src/exponential/exponential.jl` is the smallest example).
- Docstrings: signature line with `-> ReturnType`, imperative description, units in square brackets (`[kg / m³]`, `[m]`, `[rad]`), keywords under `# Keywords` with `(**Default**: ...)`. Comments are full English sentences ending in a period.
- Naming: `snake_case` functions, `_` prefix for private symbols, CamelCase types, `SCREAMING_SNAKE` constants. Keyword arguments are separated from positional ones with `;` in definitions and calls. Use explicit `return`.
- Commit messages: `<:gitmoji:> Capitalized imperative summary` up to 50 characters, textual emoji codes only (`:bug:`, `:sparkles:`, `:package:` refactor, `:racehorse:` performance, `:books:` docs, `:rotating_light:` tests, `:art:` formatting, `:wrench:` general, `:bookmark:` version), with `:warning:` appended for breaking changes. See `git log --oneline` for examples.

## Behavioral Constraints

- Model outputs are validated against reference implementations (Fortran, C, GMAT, Orekit). Never loosen a tolerance or edit a fixture to make a test pass. A change in numerical results is a bug fix or a breaking change and needs a CHANGELOG entry saying so.
- Every model validates its altitude at the entry point with `ArgumentError` (`_check_altitude`). Keep guard-clause validation at the API boundary and keep the internals assumption-free.
- The output type follows the floating-point type of the inputs (`RT = float(typeof(h))` pattern). Keep the models generic over `Number` so Float32 and dual numbers work.
- Hot paths are allocation-free and type-stable; `test/performance.jl` measures allocations and runs JET and Aqua. Prefer tuples, `StaticArraysCore` types, and named-tuple constants over heap-allocated containers.
- Every model must remain differentiable with ForwardDiff, Zygote, and Mooncake. Adding mutation or non-differentiable operations (for example `Dates` arithmetic, see `_get_doy`) requires a rule in `ext/` and coverage in `test/extensions.jl`.
- New tests follow the `@testset "Name" verbose = true begin include("./file.jl") end` wiring in `test/runtests.jl`, one file per model, with a `## Description` header explaining where the reference values came from.
- `CHANGELOG.md` gets an entry for every user-facing change: newest version first, one list item per change with a category badge (Breaking, Deprecation, Feature, Enhancement, Bugfix, Info, in that order), imperative description ending in a period, wrapped at 92 columns. Badge and PR reference links live at the end of the file.
- The docs build uses `checkdocs = :exports` and `@autodocs` in `docs/src/lib/library.md`, so every exported symbol needs a docstring. The per-model pages under `docs/src/man/` are hand-written; update them when a model's API or behavior changes.

## Not Configured

- No pre-commit hooks, no linter beyond the Aqua and JET checks inside the test suite, and no format check in CI. Do not invent one.
- No build script and no test selector, as noted above.
