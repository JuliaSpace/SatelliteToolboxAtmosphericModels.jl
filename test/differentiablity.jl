const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Enzyme", AutoEnzyme()),
    ("Mooncake", AutoMooncake(;config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
    ("Zygote", AutoZygote()),
)

##########################################################################################
# MRLMSISE-00 Fails these Packages
# 1. Enzyme: The Enzyme compiler seems to have issue with either the nested structure
# of some of the NRLMSISE-00 structures or with the type Union of the Ap space index. Either
# way it seems like it will take a bigger rework to get this working. I've tried a number of 
# attempts but end up with a large LLVM stack trace that I can't seem to resolve.
##########################################################################################
#TODO: Revisit and get this working
const _nrlmsise00_skip_backends = ["Enzyme"]

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

                @test f_fd ≈ f_ad rtol=1e-14
                @test df_fd ≈ df_ad rtol=1e-4
            end
        end
    end
end


@testset "Jacchia-Roberts 1971 Atmosphere Differentiation" begin 

    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000.0

    for backend in _BACKENDS
        if backend[1] == "Enzyme"
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
                @test f_fd ≈ f_ad rtol=1e-14
                @test df_fd ≈ something.(df_ad, 0) rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input2
                )

                # Include something() to replace Zygote "nothing" with 0.0
                @test f_fd2 ≈ f_ad2 rtol = 1e-14
                @test df_fd2 ≈ something.(df_ad2, 0) rtol=2e-1
                
            end
        end
    end
    testset_name = "JR1971 Atmosphere Enzyme"
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
                Const((x) -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density),
                AutoEnzyme(),
                input
            )
            
            # Include something() to replace Zygote "nothing" with 0.0
            @test f_fd ≈ f_ad rtol=1e-14
            @test df_fd ≈ df_ad rtol=2e-1

            f_ad2, df_ad2 = value_and_gradient(
                Const((x) -> AtmosphericModels.jr1971(x...; verbose=Val(false)).total_density),
                AutoEnzyme(),
                input2
            )

            # Include something() to replace Zygote "nothing" with 0.0
            @test f_fd2 ≈ f_ad2 rtol=1e-14
            @test df_fd2 ≈ df_ad2 rtol=2e-1
            
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

                input = [instant; h; ϕ_gd; λ; F10; F10ₐ; ap]
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


                @test f_fd ≈ f_ad rtol=1e-14
                @test df_fd ≈ df_ad atol=1e-5

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.nrlmsise00(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input2
                )

                @test f_fd2 ≈ f_ad2 rtol=1e-14
                @test df_fd2 ≈ df_ad2 atol=1e-5
            end
        end
    end
end

@testset "Jacchia-Bowman 2008 Atmosphere Differentiation" begin 

    SpaceIndices.init()

    hs = collect(90:50:1000) .* 1000.0

    for backend in _BACKENDS
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

                @test f_fd ≈ f_ad rtol=1e-14
                @test df_fd ≈ df_ad rtol=2e-1

                f_ad2, df_ad2 = value_and_gradient(
                    (x) -> AtmosphericModels.jb2008(x...; verbose=Val(false)).total_density,
                    backend[2],
                    input2
                )

                @test f_fd2 ≈ f_ad2 rtol=1e-14
                @test df_fd2 ≈ df_ad2 rtol=2e-1       
            end
        end
    end
end