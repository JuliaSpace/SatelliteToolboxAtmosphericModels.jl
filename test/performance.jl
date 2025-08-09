## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(SatelliteToolboxAtmosphericModels; ambiguities=(recursive = false), deps_compat=(check_extras = false))
end

@testset "JET Testing" begin
    rep = JET.test_package(SatelliteToolboxAtmosphericModels; toplevel_logger=nothing, target_modules=(@__MODULE__,))
end

# Verbose Val check added to logging. See link for more information.
#   https://julialang.github.io/AllocCheck.jl/dev/tutorials/optional_debugging_and_logging/

@testset "Allocation Check" begin
    rc = zeros(5)
    @test length(check_allocs(AtmosphericModels.exponential, (Float64,))) == 0

    # TODO: All of the allocations are in PolynomialRoots, this package also doesn't allow a SVector input.
    # TODO: Work through these in PolynomialRoots.jl or consider a different package.
    # TODO: 5 Allocations on Julia 1.10 CI and 10 on Julia 1 CI.

    @test length(
        check_allocs(
            (x1, x2, x3, x4, rc) -> begin
                AtmosphericModels.jr1971(
                    x1,
                    x2,
                    x3,
                    x4;
                    roots_container = rc,
                    verbose = Val(false)
                )
            end,
            (DateTime, Float64, Float64, Float64, Vector{Float64})
        )
    ) <= 10

    @test length(
        check_allocs(
            (x1, x2, x3, x4) -> begin
                AtmosphericModels.jb2008(x1, x2, x3, x4; verbose = Val(false))
            end,
            (DateTime, Float64, Float64, Float64)
        )
    ) == 0

    @test length(
        check_allocs(
            (x1, x2, x3, x4, P) -> begin
                AtmosphericModels.nrlmsise00(x1, x2, x3, x4; P = P, verbose = Val(false))
            end,
            (DateTime, Float64, Float64, Float64, Matrix{Float64})
        )
    ) == 0
end
