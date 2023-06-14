# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the Jacchia-Bowman 2008 model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] http://sol.spacenvironment.net/~JB2008/
#   [3] https://github.com/JuliaSpace/JB2008_Test
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Functions: jb2008
# ==========================================================================================

############################################################################################
#                                       Test Results
############################################################################################
#
# Scenario 01
# ===========
#
#   Compare the results with the tests in files:
#
#       JB_AUTO_OUTPUT_01.DAT
#       JB_AUTO_OUTPUT_02.DAT
#
#   which were created using JuliaSpace/JB2008_Test@ac41221. The other two files are ignored
#   due to the lack of data in the space indices.
#
############################################################################################

@testset "Fetching All Space Indices" begin
    SpaceIndices.init()

    # Scenario 01
    # ======================================================================================

    # Files with the test results.
    test_list = [
        "JB2008_AUTO_OUTPUT_01.DAT",
        "JB2008_AUTO_OUTPUT_02.DAT",
    ]

    # Execute the tests.
    for filename in test_list
        open(filename) do file
            line_num = 0
            year     = 0
            doy      = 0
            hour     = 0
            min      = 0
            sec      = 0.0

            for line in eachline(file)
                line_num += 1

                # Ignore 5 lines to skip the header.
                (line_num <= 5) && continue

                # Ignore others non important lines.
                (line_num in [7, 8, 9]) && continue

                if line_num == 6
                    # Read the next line to obtain the input data related to the
                    # time.
                    tokens = split(line)

                    year = parse(Int,     tokens[1])
                    doy  = parse(Int,     tokens[2])
                    hour = parse(Int,     tokens[3])
                    min  = parse(Int,     tokens[4])
                    sec  = parse(Float64, tokens[5])
                else
                    tokens = split(line)

                    # Read the information about the location.
                    h    = parse(Float64, tokens[1]) * 1000
                    λ    = parse(Float64, tokens[2]) |> deg2rad
                    ϕ_gd = parse(Float64, tokens[3]) |> deg2rad

                    # Read the model output.
                    exospheric_temperature = parse(Float64, tokens[4])
                    temperature = parse(Float64, tokens[5])
                    total_density = parse(Float64, tokens[6])

                    # Run the model.
                    jd  = date_to_jd(year, 1, 1, hour, min, sec) - 1 + doy
                    instant = julian2datetime(jd)
                    result = AtmosphericModels.jb2008(jd, ϕ_gd, λ, h)

                    # Compare the results.
                    @test result.exospheric_temperature ≈ exospheric_temperature atol = 0.6 rtol = 0.0
                    @test result.temperature ≈ temperature atol = 0.6 rtol = 0.0
                    @test result.total_density ≈ total_density atol = 0.0 rtol = 5e-3
                end
            end
        end
    end
end

@testset "Show" begin
    result = AtmosphericModels.jb2008(
        DateTime("2023-01-01T10:00:00"),
        0,
        0,
        500e3,
        100,
        100,
        100,
        100,
        100,
        100,
        100,
        100,
        85
    )

    expected = "JB2008 output (ρ = 3.52089e-13 kg / m³)"
    str = sprint(show, result)

    expected = """
        Jacchia-Bowman 2008 Atmospheric Model Result:
              Total density :    3.52089e-13  kg / m³
                Temperature :         923.19  K
           Exospheric Temp. :         924.79  K
          N₂ number density :    1.38441e+11  1 / m³
          O₂ number density :    3.20415e+09  1 / m³
          O  number density :    1.23957e+13  1 / m³
          Ar number density :    2.20193e+06  1 / m³
          He number density :     2.4238e+12  1 / m³
          H  number density :    4.13032e+10  1 / m³"""

    str = sprint(show, MIME("text/plain"), result)
    @test str == expected
end

@testset "Errors" begin
    @test_throws ArgumentError AtmosphericModels.jb2008(
        now(),
        0,
        0,
        89.9e3,
        100,
        100,
        100,
        100,
        100,
        100,
        100,
        100,
        85
    )
end
