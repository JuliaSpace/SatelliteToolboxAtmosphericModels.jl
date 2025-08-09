## Description #############################################################################
#
# Tests related to the Jacchia-Robert 1971 model.
#
############################################################################################

# == Functions: jr1971 =====================================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# == Scenario 01 ===========================================================================
#
#   Values obtained from GMAT R2018a using the following inputs:
#
#       Date:      2017-01-01 00:00:00 UTC
#       Latitude:  45 deg
#       Longitude: 0 deg
#       F10.7:     100
#       F10.7ₐ:    100
#       Kp:        4
#
#   Result:
#
#       | Altitude [km]  | Density [g/cm³] |
#       |----------------|------------------|
#       |             92 | 3.87506e-09      |
#       |            100 | 3.96585e-09      |
#       |          100.1 | 6.8354e-10       |
#       |          110.5 | 1.2124e-10       |
#       |            125 | 1.60849e-11      |
#       |          125.1 | 1.58997e-11      |
#       |            300 | 1.30609e-14      |
#       |            700 | 1.34785e-17      |
#       |           1500 | 4.00464e-19      |
#
# == Scenario 02 ===========================================================================
#
#   Values obtained from GMAT R2018a using the following inputs:
#
#       Date:      2017-01-01 00:00:00 UTC
#       Latitude:  45 deg
#       Longitude: 0 deg
#       F10.7:     100
#       F10.7ₐ:    100
#       Kp:        1
#
#   Result:
#
#       | Altitude [km]  | Density [g/cm³] |
#       |----------------|------------------|
#       |             92 | 3.56187e-09      |
#       |            100 | 3.65299e-09      |
#       |          100.1 | 6.28941e-10      |
#       |          110.5 | 1.11456e-10      |
#       |            125 | 1.46126e-11      |
#       |          125.1 | 1.44428e-11      |
#       |            300 | 9.08858e-15      |
#       |            700 | 8.51674e-18      |
#       |           1500 | 2.86915e-19      |
#
# == Scenario 03 ===========================================================================
#
#   Values obtained from GMAT R2018a using the following inputs:
#
#       Date:      2017-01-01 00:00:00 UTC
#       Latitude:  45 deg
#       Longitude: 0 deg
#       F10.7:     100
#       F10.7ₐ:    100
#       Kp:        9
#
#   Result:
#
#       | Altitude [km]  | Density [g/cm³] |
#       |----------------|------------------|
#       |             92 | 5.55597e-09      |
#       |            100 | 5.63386e-09      |
#       |          100.1 | 9.75634e-10      |
#       |          110.5 | 1.73699e-10      |
#       |            125 | 2.41828e-11      |
#       |          125.1 | 2.39097e-11      |
#       |            300 | 3.52129e-14      |
#       |            700 | 1.28622e-16      |
#       |           1500 | 1.97775e-18      |
#
############################################################################################

@testset "Providing All Space Indices" begin
    # Common inputs to all scenarios.
    jd      = date_to_jd(2017,1,1,0,0,0)
    instant = julian2datetime(jd)
    ϕ_gd    = 45 |> deg2rad
    λ       = 0.0
    F10     = 100.0
    F10ₐ    = 100.0
    h       = [92, 100, 100.1, 110.5, 125, 125.1, 300, 700, 1500] * 1000

    # == Scenario 01 =======================================================================

    Kp = 4

    # Results in [kg/m³].
    results = [
        3.87506e-09
        3.96585e-09
        6.83540e-10
        1.21240e-10
        1.60849e-11
        1.58997e-11
        1.30609e-14
        1.34785e-17
        4.00464e-19
    ] * 1000

    roots_container = zeros(5)

    for i in 1:length(h)
        ret = AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[i - 1 + begin], F10, F10ₐ, Kp)
        @test ret.total_density ≈ results[i - 1 + begin] rtol = 5e-4
        ret_prealloc = AtmosphericModels.jr1971(
            instant,
            ϕ_gd,
            λ,
            h[i - 1 + begin],
            F10,
            F10ₐ,
            Kp;
            roots_container = roots_container
        )
        @test ret_prealloc.total_density ≈ results[i - 1 + begin] rtol = 5e-4
    end

    # == Scenario 02 =======================================================================

    Kp = 1

    # Results in [kg/m³].
    results = [
        3.56187e-09
        3.65299e-09
        6.28941e-10
        1.11456e-10
        1.46126e-11
        1.44428e-11
        9.08858e-15
        8.51674e-18
        2.86915e-19
    ] * 1000

    for i in 1:length(h)
        ret = AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[i - 1 + begin], F10, F10ₐ, Kp)
        @test ret.total_density ≈ results[i - 1 + begin] rtol = 5e-4
        ret_prealloc = AtmosphericModels.jr1971(
            instant,
            ϕ_gd,
            λ,
            h[i - 1 + begin],
            F10,
            F10ₐ,
            Kp;
            roots_container = roots_container
        )
        @test ret_prealloc.total_density ≈ results[i - 1 + begin] rtol = 5e-4
    end

    # == Scenario 03 =======================================================================

    Kp = 9

    # Results in [kg/m³].
    results = [
        5.55597e-09
        5.63386e-09
        9.75634e-10
        1.73699e-10
        2.41828e-11
        2.39097e-11
        3.52129e-14
        1.28622e-16
        1.97775e-18
    ] * 1000

    for i in 1:length(h)
        ret = AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[i - 1 + begin], F10, F10ₐ, Kp)
        @test ret.total_density ≈ results[i - 1 + begin] rtol = 5e-4
        ret_prealloc = AtmosphericModels.jr1971(
            instant,
            ϕ_gd,
            λ,
            h[i - 1 + begin],
            F10,
            F10ₐ,
            Kp;
            roots_container = roots_container
        )
        @test ret_prealloc.total_density ≈ results[i - 1 + begin] rtol = 5e-4
    end
end

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# In this case, we already tested the function `AtmosphericModel.jr1971`. Hence, we will
# select a day and run this function with and without passing the space indices. The result
# must be the same.
#
# We have the following space indices for the instant 2023-01-01T10:00:00.000:
#
#   F10  = 152.6 sfu
#   F10ₐ = 159.12345679012347 sfu
#   Kp   = 2.667
#
############################################################################################

@testset "Fetching All Space Indices" begin
    SpaceIndices.init()

    # Expected result.
    instant = DateTime("2023-01-01T10:00:00")
    h       = collect(90:50:1000) .* 1000
    ϕ_gd    = -23 |> deg2rad
    λ       = -45 |> deg2rad
    F10     = 152.6
    F10ₐ    = 159.12345679012347
    Kp      = 2.667

    expected = AtmosphericModels.jr1971.(instant, ϕ_gd, λ, h, F10, F10ₐ, Kp)

    for k in 1:length(h)
        result = AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[k - 1 + begin])

        @test result.total_density          ≈ expected[k - 1 + begin].total_density
        @test result.temperature            ≈ expected[k - 1 + begin].temperature
        @test result.exospheric_temperature ≈ expected[k - 1 + begin].exospheric_temperature
        @test result.N2_number_density      ≈ expected[k - 1 + begin].N2_number_density
        @test result.O2_number_density      ≈ expected[k - 1 + begin].O2_number_density
        @test result.O_number_density       ≈ expected[k - 1 + begin].O_number_density
        @test result.Ar_number_density      ≈ expected[k - 1 + begin].Ar_number_density
        @test result.He_number_density      ≈ expected[k - 1 + begin].He_number_density
        @test result.H_number_density       ≈ expected[k - 1 + begin].H_number_density
    end
end

@testset "Show" begin
    result = AtmosphericModels.jr1971(
        DateTime("2023-01-01T10:00:00"),
        0,
        0,
        500e3,
        100,
        100,
        3
    )

    expected = "JR1971 output (ρ = 5.51927e-14 kg / m³)"
    str = sprint(show, result)

    expected = """
        Jacchia-Roberts 1971 Atmospheric Model Result:
              Total density :    5.15927e-14  kg / m³
                Temperature :         679.39  K
           Exospheric Temp. :         679.44  K
          N₂ number density :    1.47999e+09  1 / m³
          O₂ number density :    1.50981e+07  1 / m³
          O  number density :    1.58322e+12  1 / m³
          Ar number density :        2149.38  1 / m³
          He number density :    1.42329e+12  1 / m³
          H  number density :              0  1 / m³"""

    str = sprint(show, MIME("text/plain"), result)
    @test str == expected
end

@testset "Errors" begin
    @test_throws ArgumentError AtmosphericModels.jr1971(now(), 0, 0, 89.9e3, 100, 100, 3)
    roots_container = zeros(6)
    @test_throws ArgumentError AtmosphericModels.jr1971(
        now(),
        0,
        0,
        89.9e3,
        100,
        100,
        3;
        roots_container = roots_container
    )
end
