## Description #############################################################################
#
# Tests to verify that package extensions load correctly and provide expected functionality.
#
############################################################################################

@testset "ForwardDiff Differentiation" begin
    SpaceIndices.init()

    instant = datetime2julian(DateTime("2023-01-01T10:00:00"))
    ϕ_gd    = deg2rad(-23.0)
    λ       = deg2rad(-45.0)

    # We test altitudes in every branch of the model, including the region below 125 km in
    # which the roots of the quartic polynomial are required.
    for h in (95e3, 110e3, 500e3)
        g = ForwardDiff.gradient(
            x -> AtmosphericModels.jr1971(x...; verbose = Val(false)).total_density,
            [instant, ϕ_gd, λ, h],
        )
        @test all(isfinite, g)
        @test !all(iszero, g)
    end
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

    # Compute the gradient using Zygote and compare it against ForwardDiff. This exercises
    # the pullback, which must return one tangent per positional argument.
    g_zygote = Zygote.gradient(
        (jd, h, ϕ_gd, λ, F10ₐ, F10, ap) ->
            AtmosphericModels.nrlmsise00(jd, h, ϕ_gd, λ, F10ₐ, F10, ap).total_density,
        jd,
        h,
        ϕ_gd,
        λ,
        F10ₐ,
        F10,
        ap,
    )

    g_forwarddiff = ForwardDiff.gradient(
        x -> AtmosphericModels.nrlmsise00(x...).total_density,
        [jd, h, ϕ_gd, λ, F10ₐ, F10, ap],
    )

    for i in 1:7
        @test g_zygote[i] ≈ g_forwarddiff[i] rtol = 1e-10
    end

    # Differentiating a call with a pre-allocated Legendre matrix must also work since the
    # pullback must not forward the Float64 matrix to the dual-valued evaluation.
    P = zeros(9, 9)

    g_prealloc = Zygote.gradient(
        (jd, h, ϕ_gd, λ, F10ₐ, F10, ap) ->
            AtmosphericModels.nrlmsise00(jd, h, ϕ_gd, λ, F10ₐ, F10, ap; P = P).total_density,
        jd,
        h,
        ϕ_gd,
        λ,
        F10ₐ,
        F10,
        ap,
    )

    for i in 1:7
        @test g_prealloc[i] ≈ g_forwarddiff[i] rtol = 1e-10
    end
end
