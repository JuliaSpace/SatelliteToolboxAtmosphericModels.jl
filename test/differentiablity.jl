## Description #############################################################################
#
# Tests related to the automatic differentiation of the atmospheric models.
#
############################################################################################

const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
)

@testset "Exponential Atmosphere Differentiation" begin
    hs = collect(90:50:1000) .* 1000.0
    for h in hs
        f_fd, df_fd = value_and_derivative(
            (x) -> AtmosphericModels.exponential(x),
            AutoFiniteDiff(),
            h
        )

        for backend in _BACKENDS
            @eval @testset $("Exponential Atmosphere h=$h m " * string(backend[1])) begin
                f_ad, df_ad = value_and_derivative(
                    (x) -> AtmosphericModels.exponential(x),
                    $backend[2],
                    $h
                )

                @test $f_fd == f_ad
                @test $df_fd ≈ df_ad rtol=1e-4
            end
        end
    end
end


@testset "Jacchia-Roberts 1971 Atmosphere Differentiation" begin
    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000
    for h in hs
        instant = datetime2julian(DateTime("2023-01-01T10:00:00"))
        ϕ_gd    = -23 |> deg2rad
        λ       = -45 |> deg2rad
        F10     = 152.6
        F10ₐ    = 159.12345679012347
        Kp      = 2.667

        input = [instant; ϕ_gd; λ; h; F10; F10ₐ; Kp]

        f_fd, df_fd = value_and_gradient(
            (x) -> AtmosphericModels.jr1971(x...).total_density,
            AutoFiniteDiff(),
            input
        )

        f_fd2, df_fd2 = value_and_gradient(
            (x) -> AtmosphericModels.jr1971(x...).total_density,
            AutoFiniteDiff(),
            input[1:4]
        )

        for backend in _BACKENDS
            @eval @testset $("JR1971 Atmosphere h=$h m " * string(backend[1])) begin
                f_ad, df_ad = value_and_gradient(
                    (x) -> AtmosphericModels.jr1971(x...).total_density,
                    $backend[2],
                    $input
                )

                @test $f_fd == f_ad
                @test $df_fd ≈ df_ad rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.jr1971(x...).total_density,
                    $backend[2],
                    $input[1:4]
                )

                @test $f_fd2 == f_ad2
                @test $df_fd2 ≈ df_ad2 rtol=2e-1
            end
        end
    end
end


@testset "NRLMSISE-00 Atmosphere Differentiation" begin
    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000
    for h in hs
        instant = datetime2julian(DateTime("2023-01-01T10:00:00"))
        ϕ_gd    = -23 |> deg2rad
        λ       = -45 |> deg2rad
        F10     = 121
        F10ₐ    = 80
        ap      = 7

        input = [instant; ϕ_gd; λ; h; F10; F10ₐ; ap]

        f_fd, df_fd = value_and_gradient(
            (x) -> AtmosphericModels.nrlmsise00(x...).total_density,
            AutoFiniteDiff(),
            input
        )

        f_fd2, df_fd2 = value_and_gradient(
            (x) -> AtmosphericModels.nrlmsise00(x...).total_density,
            AutoFiniteDiff(),
            input[1:4]
        )

        for backend in _BACKENDS
            @eval @testset $("NRLMSISE-00 Atmosphere h=$h m" * string(backend[1])) begin
                f_ad, df_ad = value_and_gradient(
                    (x) -> AtmosphericModels.nrlmsise00(x...).total_density,
                    $backend[2],
                    $input
                )

                @test $f_fd == f_ad
                @test $df_fd ≈ df_ad rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.nrlmsise00(x...).total_density,
                    $backend[2],
                    $input[1:4]
                )

                @test $f_fd2 == f_ad2
                @test $df_fd2 ≈ df_ad2 rtol=2e-1
            end
        end
    end
end

@testset "Jacchia-Bowman 2008 Atmosphere Differentiation" begin
    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000

    for h in hs
        input = [datetime2julian(DateTime("2023-01-01T10:00:00")); 0; 0; 500e3; 100; 100; 100; 100; 100; 100; 100; 100; 85]

        f_fd, df_fd = value_and_gradient(
            (x) -> AtmosphericModels.jb2008(x...).total_density,
            AutoFiniteDiff(),
            input
        )

        f_fd2, df_fd2 = value_and_gradient(
            (x) -> AtmosphericModels.jb2008(x...).total_density,
            AutoFiniteDiff(),
            input[1:4]
        )

        for backend in _BACKENDS
            @eval @testset $("JB2008 Atmosphere h=$h m" * string(backend[1])) begin
                f_ad, df_ad = value_and_gradient(
                    (x) -> AtmosphericModels.jb2008(x...).total_density,
                    $backend[2],
                    $input
                )

                @test $f_fd == f_ad
                @test $df_fd ≈ df_ad rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.jb2008(x...).total_density,
                    $backend[2],
                    $input[1:4]
                )

                @test $f_fd2 == f_ad2
                @test $df_fd2 ≈ df_ad2 rtol=2e-1
            end
        end
    end
end
