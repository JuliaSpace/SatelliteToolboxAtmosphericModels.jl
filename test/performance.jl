## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(SatelliteToolboxAtmosphericModels; ambiguities = (recursive = false))
end

if VERSION >= v"1.12"
    @warn "JET.jl test skipped on Julia 1.12+ due to MethodTableView incompatibility"
else
    @testset "JET Testing" begin
        rep = JET.test_package(
            SatelliteToolboxAtmosphericModels;
            toplevel_logger = nothing,
            target_modules = (SatelliteToolboxAtmosphericModels,),
        )
    end
end

############################################################################################
#                                   Runtime Allocations                                    #
############################################################################################

# Measure the allocations of a call after compiling it. The function and its arguments
# are passed explicitly to avoid capturing global variables in closures, which allocate.
# The number of arguments is a type parameter so that the splatted call is inferred on
# Julia 1.10, where an untyped `args...` boxes the returned value, reporting spurious
# allocations.
function _measure_allocations(f::F, args::Vararg{Any, N}) where {F, N}
    f(args...)
    return @allocated f(args...)
end

@testset "Runtime Allocations" begin
    SpaceIndices.init()

    jd   = datetime2julian(DateTime("2023-01-01T10:00:00"))
    ϕ_gd = deg2rad(-23.0)
    λ    = deg2rad(-45.0)
    P    = zeros(8, 4)

    # == Methods With Explicit Space Indices ===============================================
    #
    # These methods are used in the hot paths of the orbit propagators and must not
    # allocate.

    @test _measure_allocations(AtmosphericModels.exponential, 300e3) == 0
    @test _measure_allocations(AtmosphericModels.harrispriester, jd, ϕ_gd, λ, 300e3) == 0
    @test _measure_allocations(
        AtmosphericModels.harrispriester_modified, jd, ϕ_gd, λ, 300e3, 150.0
    ) == 0

    for h in (95e3, 110e3, 300e3, 700e3)
        @test _measure_allocations(
            AtmosphericModels.jr1971, jd, ϕ_gd, λ, h, 120.0, 118.0, 3.0
        ) == 0
        @test _measure_allocations(
            AtmosphericModels.jb2008,
            jd,
            ϕ_gd,
            λ,
            h,
            120.0,
            118.0,
            100.0,
            99.0,
            100.0,
            99.0,
            100.0,
            99.0,
            20.0,
        ) == 0
        @test _measure_allocations(
            AtmosphericModels.jacchia1977, jd, ϕ_gd, λ, h, 120.0, 118.0, 3.0
        ) == 0
        @test _measure_allocations(
            (a1, a2, a3, a4, a5, a6, a7) -> AtmosphericModels.jacchia1977(
                a1, a2, a3, a4, a5, a6, a7; variant = Val(:stela)
            ),
            jd,
            ϕ_gd,
            λ,
            h,
            120.0,
            118.0,
            3.0,
        ) == 0
    end

    for h in (50e3, 100e3, 300e3)
        @test _measure_allocations(
            (a1, a2, a3, a4, a5, a6, a7) ->
                AtmosphericModels.nrlmsise00(a1, a2, a3, a4, a5, a6, a7; P = P),
            jd,
            h,
            ϕ_gd,
            λ,
            118.0,
            120.0,
            10.0,
        ) == 0
        @test _measure_allocations(
            (a1, a2, a3, a4, a5, a6, a7) ->
                AtmosphericModels.nrlmsise00(a1, a2, a3, a4, a5, a6, a7; P = P),
            jd,
            h,
            ϕ_gd,
            λ,
            118.0,
            120.0,
            [10.0, 12.0, 9.0, 11.0, 10.0, 8.0, 7.0],
        ) == 0
    end

    # == Methods With Automatic Space Index Fetching =======================================
    #
    # These methods must not allocate besides the debug message, which is only built if
    # the debug logging is enabled.

    @test _measure_allocations(AtmosphericModels.jr1971, jd, ϕ_gd, λ, 300e3) == 0
    @test _measure_allocations(AtmosphericModels.jb2008, jd, ϕ_gd, λ, 300e3) == 0
    @test _measure_allocations(AtmosphericModels.jacchia1977, jd, ϕ_gd, λ, 300e3) == 0
    @test _measure_allocations(
        AtmosphericModels.harrispriester_modified, jd, ϕ_gd, λ, 300e3
    ) == 0
    @test _measure_allocations(
        (a1, a2, a3, a4) -> AtmosphericModels.nrlmsise00(a1, a2, a3, a4; P = P),
        jd,
        300e3,
        ϕ_gd,
        λ,
    ) == 0
end

############################################################################################
#                                    Static Allocations                                    #
############################################################################################

# AllocCheck.jl proves statically that the compiled code of the methods with explicit space
# indices cannot allocate. The check is skipped on Julia 1.12+ because AllocCheck.jl
# detects runtime calls of the newer Julia versions (e.g. jl_get_pgcstack_static) as
# allocations. The runtime allocation tests above cover all the Julia versions. The
# methods with automatic space index fetching are not checked here since the debug message
# they emit is reported as a potential allocation even when the logging is disabled.
if VERSION >= v"1.12"
    @warn "Static allocation tests skipped on Julia 1.12+ (AllocCheck.jl limitation)."
else
    @testset "Allocation Check" begin
        @test length(check_allocs(AtmosphericModels.exponential, (Float64,))) == 0

        @test length(
            check_allocs(
                AtmosphericModels.jr1971,
                (DateTime, Float64, Float64, Float64, Float64, Float64, Float64),
            ),
        ) == 0

        @test length(
            check_allocs(
                AtmosphericModels.jacchia1977,
                (DateTime, Float64, Float64, Float64, Float64, Float64, Float64),
            ),
        ) == 0

        @test length(
            check_allocs(
                AtmosphericModels.jb2008,
                (
                    DateTime,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                ),
            ),
        ) == 0

        @test length(
            check_allocs(
                (x1, x2, x3, x4, x5, x6, x7, P) -> begin
                    AtmosphericModels.nrlmsise00(x1, x2, x3, x4, x5, x6, x7; P = P)
                end,
                (
                    DateTime,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Matrix{Float64},
                ),
            ),
        ) == 0

        @test length(
            check_allocs(
                AtmosphericModels.harrispriester, (DateTime, Float64, Float64, Float64)
            ),
        ) == 0

        @test length(
            check_allocs(
                AtmosphericModels.harrispriester_modified,
                (DateTime, Float64, Float64, Float64, Float64),
            ),
        ) == 0
    end
end
