## Description #############################################################################
#
# Tests to verify that package extensions load correctly and provide expected functionality.
#
############################################################################################

@testset "ForwardDiff Extension" begin
    SpaceIndices.init()

    instant = datetime2julian(DateTime("2023-01-01T10:00:00"))
    ϕ_gd    = deg2rad(-23.0)
    λ       = deg2rad(-45.0)
    h       = 500e3

    g = ForwardDiff.gradient(
        x -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density,
        [instant, ϕ_gd, λ, h]
    )
    @test all(isfinite, g)
    @test !all(iszero, g)
end

@testset "Mooncake Extension" begin
    jd = datetime2julian(DateTime("2023-01-01T10:00:00"))

    y, pb = ChainRulesCore.rrule(AtmosphericModels._get_doy, jd)
    @test y == AtmosphericModels._get_doy(jd)
    nt, dy = pb(1.0)
    @test nt isa ChainRulesCore.NoTangent
    @test dy == 1.0
end

@testset "Zygote Extension" begin
    SpaceIndices.init()

    jd   = datetime2julian(DateTime("2023-01-01T10:00:00"))
    h    = 500e3
    ϕ_gd = deg2rad(-23.0)
    λ    = deg2rad(-45.0)
    F10ₐ = 80.0
    F10  = 121.0
    ap   = 7.0

    result, pullback = Zygote.ChainRulesCore.rrule(
        AtmosphericModels.nrlmsise00, jd, h, ϕ_gd, λ, F10ₐ, F10, ap
    )
    @test result.total_density > 0
end
