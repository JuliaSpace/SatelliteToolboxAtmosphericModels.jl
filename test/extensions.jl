## Description #############################################################################
#
# Tests to verify that the package extensions load correctly and provide expected
# functionality with the automatic differentiation backends.
#
############################################################################################

@testset "ForwardDiff Differentiation" begin
    SpaceIndices.init()

    jd   = datetime2julian(DateTime("2023-01-01T10:00:00"))
    ϕ_gd = deg2rad(-23.0)
    λ    = deg2rad(-45.0)

    # We test altitudes in every branch of the JR1971 model, including the region below
    # 100 km integrated numerically and the region below 125 km in which the roots of the
    # quartic polynomial are required.
    for h in (95e3, 110e3, 500e3)
        g = ForwardDiff.gradient(
            x -> AtmosphericModels.jr1971(x...).total_density, [jd, ϕ_gd, λ, h]
        )
        @test all(isfinite, g)
        @test !all(iszero, g)
    end

    # All the other models must be differentiable with respect to the position and epoch
    # using the explicit-index methods.
    models = (
        x -> AtmosphericModels.jb2008(
            x[1], x[2], x[3], x[4], 120, 118, 100, 99, 100, 99, 100, 99, 20
        ).total_density,
        x -> AtmosphericModels.jacchia1977(
            x[1], x[2], x[3], x[4], 120, 118, 3
        ).total_density,
        x -> AtmosphericModels.jacchia1977(
            x[1], x[2], x[3], x[4], 120, 118, 3; variant = Val(:stela)
        ).total_density,
        x -> AtmosphericModels.nrlmsise00(
            x[1], x[4], x[2], x[3], 118, 120, 10
        ).total_density,
        x -> AtmosphericModels.harrispriester(x[1], x[2], x[3], x[4]),
        x -> AtmosphericModels.harrispriester_modified(x[1], x[2], x[3], x[4], 150),
    )

    for f in models, h in (150e3, 300e3)
        g = ForwardDiff.gradient(f, [jd, ϕ_gd, λ, h])
        @test all(isfinite, g)
        @test !all(iszero, g)
    end

    # The exponential model depends only on the altitude.
    for h in (0.0, 300e3, 1500e3)
        d = ForwardDiff.derivative(AtmosphericModels.exponential, h)
        @test isfinite(d)
        @test d < 0
    end
end

@testset "ChainRulesCore Extension" begin
    jd = datetime2julian(DateTime("2023-01-01T10:00:00"))

    y, pb = ChainRulesCore.rrule(AtmosphericModels._get_doy, jd)
    @test y == AtmosphericModels._get_doy(jd)
    nt, dy = pb(1.0)
    @test nt isa ChainRulesCore.NoTangent
    @test dy == 1.0
end

@testset "Mooncake Extension" begin
    jd   = datetime2julian(DateTime("2023-01-01T10:00:00"))
    ϕ_gd = deg2rad(-23.0)
    λ    = deg2rad(-45.0)

    # The reverse-mode gradients must match the forward-mode ones. The day of the year is
    # computed with `Dates` arithmetic, whose derivative is provided by the extension.
    models = (
        jd -> AtmosphericModels.jr1971(jd, ϕ_gd, λ, 300e3, 120.0, 118.0, 3.0).total_density,
        jd -> AtmosphericModels.nrlmsise00(
            jd, 300e3, ϕ_gd, λ, 118.0, 120.0, 10.0
        ).total_density,
    )

    for f in models
        rule    = Mooncake.build_rrule(f, jd)
        v, grad = Mooncake.value_and_gradient!!(rule, f, jd)

        @test v == f(jd)
        @test grad[2] ≈ ForwardDiff.derivative(f, jd) rtol = 1e-8
    end
end

@testset "Zygote Extension" begin
    jd   = datetime2julian(DateTime("2023-01-01T10:00:00"))
    h    = 500e3
    ϕ_gd = deg2rad(-23.0)
    λ    = deg2rad(-45.0)
    F10ₐ = 80.0
    F10  = 121.0

    # == Daily Magnetic Index ==============================================================

    ap = 7.0

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
    P = zeros(8, 4)

    g_prealloc = Zygote.gradient(
        (jd, h, ϕ_gd, λ, F10ₐ, F10, ap) ->
            AtmosphericModels.nrlmsise00(
                jd, h, ϕ_gd, λ, F10ₐ, F10, ap; P = P
            ).total_density,
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

    # == Magnetic Index Vector =============================================================

    ap = [7.0, 9.0, 12.0, 6.0, 5.0, 8.0, 7.0]

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
        x -> AtmosphericModels.nrlmsise00(
            x[1], x[2], x[3], x[4], x[5], x[6], x[7:13]
        ).total_density,
        vcat([jd, h, ϕ_gd, λ, F10ₐ, F10], ap),
    )

    for i in 1:6
        @test g_zygote[i] ≈ g_forwarddiff[i] rtol = 1e-10
    end

    @test g_zygote[7] isa AbstractVector
    @test length(g_zygote[7]) == 7

    for i in 1:7
        @test g_zygote[7][i] ≈ g_forwarddiff[6 + i] rtol = 1e-10
    end
end
