const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Diffractor", AutoDiffractor()),
    ("Enzyme", AutoEnzyme()),
    ("Mooncake", AutoMooncake(;config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
    ("Zygote", AutoZygote()),
)

##########################################################################################
# JR1971 Fails these Packages
# 1. Diffractor: Doesn't Like how the Average Space Index is Calculated
    # This is likely on Diffractors end, an issue will be raised and these tests can be updated once resolved
# 2. Enzyme: Doesn't Like the Polynomial roots computation, could potentially look at other roots packages, writing a custom rule
    # This is hard to debug as it causes a segfault
##########################################################################################
const _jr1971_skip_backends = ["Diffractor", "Enzyme"]

##########################################################################################
# JB2008 Fails these Packages
# 1. Diffractor: Doesn't Like how the Average Space Index is Calculated
    # This is likely on Diffractors end, an issue will be raised and these tests can be updated once resolved
    # ALso an error in the bitcast of the _jb2008_∫, could look into replacing this with Integrals.jl
##########################################################################################
const _jb2008_skip_backends = ["Diffractor"]

##########################################################################################
# MRLMSISE-00 Fails these Packages
# 1. Diffractor: Doesn't Like how the Average Space Index is Calculated
    # This is likely on Diffractors end, an issue will be raised and these tests can be updated once resolved
# 2. Zygote: The fail due to the legendre computation being mutating.
    # This change would need to happen in SatelliteToolboxLegendre.jl
# 3. Enzyme: Stack Overflow Error, I don't know why this is happening
##########################################################################################
const _nrlmsise00_skip_backends = ["Diffractor", "Enzyme", "Zygote"]

@testset "Exponential Atmosphere Differentiation" begin
    
    hs = collect(90:50:1000) .* 1000.0

    for backend in _BACKENDS
        testset_name = "Exponential Atmosphere " * string(backend[1])
        @testset "$testset_name" begin
            for h in hs
                f_fd, df_fd = value_and_derivative(
                    (x) -> AtmosphericModels.exponential(x),
                    AutoFiniteDiff(),
                    h
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> AtmosphericModels.exponential(x),
                    backend[2],
                    h
                )

                @test f_fd == f_ad
                @test df_fd ≈ df_ad rtol=1e-4
            end
        end
    end
end


@testset "Jacchia-Roberts 1971 Atmosphere Differentiation" begin 

    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000.0

    for backend in _BACKENDS
        if backend[1] in _jr1971_skip_backends
            continue
        end
        testset_name = "JR1971 Atmosphere " * string(backend[1])
        @testset "$testset_name" begin
            for h in hs
                instant = datetime2julian(DateTime("2023-01-01T10:00:00"))
                ϕ_gd    = -23 |> deg2rad
                λ       = -45 |> deg2rad
                F10     = 152.6
                F10ₐ    = 159.12345679012347
                Kp      = 2.667

                input = [instant; ϕ_gd; λ; h; F10; F10ₐ; Kp]
                input2 = input[1:4]

                f_fd, df_fd = value_and_gradient(
                    (x) -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density,
                    AutoFiniteDiff(),
                    input
                )

                f_fd2, df_fd2 = value_and_gradient(
                    (x) -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density,
                    AutoFiniteDiff(),
                    input2
                )

                
                f_ad, df_ad = value_and_gradient(
                    (x) -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input
                )
                
                # Include something() to replace Zygote "nothing" with 0.0
                @test f_fd == f_ad
                @test df_fd ≈ something.(df_ad, 0) rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input2
                )

                # Include something() to replace Zygote "nothing" with 0.0
                @test f_fd2 == f_ad2
                @test df_fd2 ≈ something.(df_ad2, 0) rtol=2e-1
                
            end
        end
    end
end


@testset "NRLMSISE-00 Atmosphere Differentiation" begin 

    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000.0

    for backend in _BACKENDS
        if backend[1] in _nrlmsise00_skip_backends
            continue
        end
        testset_name = "NRLMSISE-00 Atmosphere " * string(backend[1])
        @testset "$testset_name" begin
            for h in hs
                instant = datetime2julian(DateTime("2023-01-01T10:00:00"))
                ϕ_gd    = -23 |> deg2rad
                λ       = -45 |> deg2rad
                F10     = 121
                F10ₐ    = 80
                ap      = 7

                input = [instant; ϕ_gd; λ; h; F10; F10ₐ; ap]
                input2 = input[1:4]


                f_fd, df_fd = value_and_gradient(
                    (x) -> AtmosphericModels.nrlmsise00(x...).total_density,
                    AutoFiniteDiff(),
                    input
                )

                f_fd2, df_fd2 = value_and_gradient(
                    (x) -> AtmosphericModels.nrlmsise00(x...; verbose=Val(false)).total_density,
                    AutoFiniteDiff(),
                    input2
                )

                f_ad, df_ad = value_and_gradient(
                    (x) -> AtmosphericModels.nrlmsise00(x...).total_density,
                    backend[2],
                    input
                )

                @test f_fd == f_ad
                @test df_fd ≈ df_ad rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.nrlmsise00(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input2
                )

                @test f_fd2 == f_ad2
                @test df_fd2 ≈ df_ad2 rtol=2e-1
            end
        end
    end
end

@testset "Jacchia-Bowman 2008 Atmosphere Differentiation" begin 

    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000.0

    for backend in _BACKENDS
        if backend[1] in _jb2008_skip_backends
            continue
        end
        testset_name = "JB2008 Atmosphere " * string(backend[1])
        @testset "$testset_name" begin
            for h in hs
                input = [datetime2julian(DateTime("2023-01-01T10:00:00")); 0; 0; 500e3; 100; 100; 100; 100; 100; 100; 100; 100; 85]
                input2 = input[1:4]

                f_fd, df_fd = value_and_gradient(
                    (x) -> AtmosphericModels.jb2008(x...; verbose=Val(true)).total_density,
                    AutoFiniteDiff(),
                    input
                )

                f_fd2, df_fd2 = value_and_gradient(
                    (x) -> AtmosphericModels.jb2008(x...; verbose=Val(true)).total_density,
                    AutoFiniteDiff(),
                    input2
                )

                f_ad, df_ad = value_and_gradient(
                    (x) -> AtmosphericModels.jb2008(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input
                )

                @test f_fd == f_ad
                @test df_fd ≈ df_ad rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.jb2008(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input2
                )

                @test f_fd2 == f_ad2
                @test df_fd2 ≈ df_ad2 rtol=2e-1       
            end
        end
    end
end