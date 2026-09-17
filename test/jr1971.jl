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
# NOTE: The values at 92 km and 100 km in the tables below were obtained with this
# implementation after fixing the density between 90 km and 100 km, and not from GMAT
# R2018a. GMAT reproduces the closed-form solution of the reference, which yields an almost
# constant density in this region (e.g. 3.87506e-09 g/cm³ at 92 km and 3.96585e-09 g/cm³ at
# 100 km in the scenario 01, both higher than the density at 90 km) and a discontinuity of a
# factor of about 6 at 100 km. The corrected values are validated by the numerical
# integration of the barometric equation in the testset "Density Between 90 km and 100 km".
# They are marked with (*) below.
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
#       |             92 | 2.70372e-09 (*)  |
#       |            100 | 6.95226e-10 (*)  |
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
#       |             92 | 2.48519e-09 (*)  |
#       |            100 | 6.39684e-10 (*)  |
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
#       |             92 | 3.87664e-09 (*)  |
#       |            100 | 9.92375e-10 (*)  |
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
    jd      = date_to_jd(2017, 1, 1, 0, 0, 0)
    instant = julian2datetime(jd)
    ϕ_gd    = 45 |> deg2rad
    λ       = 0.0
    F10     = 100.0
    F10ₐ    = 100.0
    h       = [92, 100, 100.1, 110.5, 125, 125.1, 300, 700, 1500] * 1000

    # == Scenario 01 =======================================================================

    Kp = 4

    # Results in [kg/m³].
    results =
        [
            2.70372e-09
            6.95226e-10
            6.83540e-10
            1.21240e-10
            1.60849e-11
            1.58997e-11
            1.30609e-14
            1.34785e-17
            4.00464e-19
        ] * 1000

    for i in 1:length(h)
        ret = AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[i - 1 + begin], F10, F10ₐ, Kp)
        @test ret.total_density ≈ results[i - 1 + begin] rtol = 5e-4
    end

    # == Scenario 02 =======================================================================

    Kp = 1

    # Results in [kg/m³].
    results =
        [
            2.48519e-09
            6.39684e-10
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
    end

    # == Scenario 03 =======================================================================

    Kp = 9

    # Results in [kg/m³].
    results =
        [
            3.87664e-09
            9.92375e-10
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
    end
end

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# The nighttime minimum global exospheric temperature is:
#
#   Tc = 379 + 3.24 F10ₐ + 1.3 (F10 - F10ₐ),
#
# where the term multiplied by 3.24 must use the 81-day averaged flux F10ₐ [1, 2]. All
# GMAT-based scenarios above use F10 == F10ₐ and, thus, cannot detect an error in this
# term. The values below are regression snapshots obtained from this implementation after
# fixing the formula, using F10 = 150, F10ₐ = 100, and Kp = 4. For reference, the buggy
# formula (3.24 F10) leads to T∞ ≈ 1056 K at 300 km, whereas the correct one yields
# T∞ ≈ 894.35 K. The snapshot at 100 km was updated after fixing the density between 90 km
# and 100 km.
#
############################################################################################

@testset "Density Between 90 km and 100 km" begin
    # The closed-form solution of the barometric equation between 90 km and 100 km in the
    # reference [1] (and in GMAT) leads to a density almost constant in this region and a
    # discontinuity of a factor of about 6 at 100 km. The model now integrates the
    # barometric equation numerically. Hence, the density must decrease monotonically and
    # must be continuous at the limits of the region.
    jd   = date_to_jd(2017, 1, 1, 0, 0, 0)
    ϕ_gd = deg2rad(45)
    λ    = 0.0

    for Kp in (0, 4, 9), F10 in (70.0, 150.0, 250.0)
        ρ = [
            AtmosphericModels.jr1971(jd, ϕ_gd, λ, h, F10, F10, Kp).total_density for
            h in 90e3:500:100e3
        ]

        @test all(diff(ρ) .< 0)

        # The density must drop by a factor between 4 and 8 in this region.
        @test 4 < ρ[1] / ρ[end] < 8

        # Continuity at 90 km and 100 km (the tolerance accounts for the density gradient
        # and for the polynomial fit of the density at 100 km used above this altitude).
        for h in (90e3, 100e3)
            ρ₋ = AtmosphericModels.jr1971(jd, ϕ_gd, λ, h, F10, F10, Kp).total_density
            ρ₊ = AtmosphericModels.jr1971(jd, ϕ_gd, λ, h + 1e-3, F10, F10, Kp).total_density
            @test ρ₊ ≈ ρ₋ rtol = 2e-3
        end
    end

    # The 8-point Gauss-Legendre quadrature of the integrand of the barometric equation
    # must match a fine midpoint integration between 90 km and 100 km.
    out = AtmosphericModels.jr1971(jd, ϕ_gd, λ, 95e3, 100.0, 100.0, 4)
    T∞  = out.exospheric_temperature
    Tx  = 371.6678 + 0.0518806 * T∞ - 294.3505 * exp(-0.00216222 * T∞)
    C   = AtmosphericModels._JR1971_CONSTANTS

    f(z) =
        C.g₀ * C.Ra^2 / (C.Ra + z)^2 * AtmosphericModels._jr1971_mean_molecular_mass(z) /
        AtmosphericModels._jr1971_temperature(z, Tx, T∞) / C.Rstar

    n  = 200_000
    Δh = 10 / n
    I  = sum(f(90 + (i + 0.5) * Δh) for i in 0:(n - 1)) * Δh

    G = 5 * sum(
        AtmosphericModels._GAUSS_LEGENDRE_8_WEIGHTS[i] *
        f(95 + 5 * AtmosphericModels._GAUSS_LEGENDRE_8_NODES[i]) for i in 1:8
    )

    @test G ≈ I rtol = 1e-9
end

@testset "Exospheric Temperature When F10 != F10ₐ" begin
    jd      = date_to_jd(2017, 1, 1, 0, 0, 0)
    instant = julian2datetime(jd)
    ϕ_gd    = 45 |> deg2rad
    λ       = 0.0
    F10     = 150.0
    F10ₐ    = 100.0
    Kp      = 4.0

    h = [100, 125.1, 300, 700, 1500] * 1000

    expected_ρ = [
        6.941830425832428e-7
        1.615563485591969e-8
        1.647659416350587e-11
        2.0027542699016525e-14
        5.588020206270033e-16
    ]

    expected_T∞ = [
        837.8040529227294
        837.8040529227294
        894.3500344230608
        894.3500344230608
        894.3500344230608
    ]

    for i in 1:length(h)
        ret = AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[i - 1 + begin], F10, F10ₐ, Kp)
        @test ret.total_density ≈ expected_ρ[i - 1 + begin] rtol = 1e-6
        @test ret.exospheric_temperature ≈ expected_T∞[i - 1 + begin] rtol = 1e-6
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
# We have the following space indices for the instant 2023-01-01T10:00:00.000, using the
# 10.7-cm flux adjusted to 1 AU, as used to fit the Jacchia models:
#
#   F10  = 147.5 sfu
#   F10ₐ = 154.2407407407407 sfu
#   Kp   = 2.0 (3-hour delayed value, i.e. related to the interval 06:00 - 09:00)
#
############################################################################################

@testset "Fetching All Space Indices" begin
    SpaceIndices.init()

    # Expected result.
    instant = DateTime("2023-01-01T10:00:00")
    h       = collect(90:50:1000) .* 1000
    ϕ_gd    = -23 |> deg2rad
    λ       = -45 |> deg2rad
    F10     = 147.5
    F10ₐ    = 154.2407407407407
    Kp      = 2.0

    expected = AtmosphericModels.jr1971.(instant, ϕ_gd, λ, h, F10, F10ₐ, Kp)

    for k in 1:length(h)
        result = AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[k - 1 + begin])

        @test result.total_density ≈ expected[k - 1 + begin].total_density
        @test result.temperature ≈ expected[k - 1 + begin].temperature
        @test result.exospheric_temperature ≈ expected[k - 1 + begin].exospheric_temperature
        @test result.N2_number_density ≈ expected[k - 1 + begin].N2_number_density
        @test result.O2_number_density ≈ expected[k - 1 + begin].O2_number_density
        @test result.O_number_density ≈ expected[k - 1 + begin].O_number_density
        @test result.Ar_number_density ≈ expected[k - 1 + begin].Ar_number_density
        @test result.He_number_density ≈ expected[k - 1 + begin].He_number_density
        @test result.H_number_density ≈ expected[k - 1 + begin].H_number_density
    end

    # == Day Boundary ======================================================================
    #
    # At 2023-01-02T01:00:00, the 3-hour delayed instant is 2023-01-01T22:00:00. Hence, the
    # Kp must be taken from the last 3-hour interval of 2023-01-01 (Kp = 4.0) instead of
    # the first interval of 2023-01-02. The daily and averaged F10.7 for this instant are
    # the same as in the previous test.

    instant_boundary = DateTime("2023-01-02T01:00:00")
    expected_boundary = AtmosphericModels.jr1971(
        instant_boundary, ϕ_gd, λ, 300e3, F10, F10ₐ, 4.0
    )
    result_boundary = AtmosphericModels.jr1971(instant_boundary, ϕ_gd, λ, 300e3)

    @test result_boundary.total_density ≈ expected_boundary.total_density
    @test result_boundary.temperature ≈ expected_boundary.temperature
end

@testset "Show" begin
    result = AtmosphericModels.jr1971(
        DateTime("2023-01-01T10:00:00"), 0, 0, 500e3, 100, 100, 3
    )

    expected = "JR1971 output (ρ = 5.51927e-14 kg / m³)"
    str = sprint(show, result)

    expected = """
          Jacchia-Roberts 1971 Atmospheric Model Result:
                Total density :    3.63066e-13  kg / m³
                  Temperature :         907.60  K
             Exospheric Temp. :         907.92  K
            N₂ number density :    7.30485e+10  1 / m³
            O₂ number density :    1.35105e+09  1 / m³
            O  number density :    1.29933e+13  1 / m³
            Ar number density :         630285  1 / m³
            He number density :    2.16503e+12  1 / m³
            H  number density :              0  1 / m³"""

    str = sprint(show, MIME("text/plain"), result)
    @test str == expected
end

@testset "Helium Seasonal Correction at the Equinoxes" begin
    # The correction must be 0 (and not NaN) when the Sun declination is exactly zero. The
    # previous formulation contained the term δs / (2 abs(δs)), which is NaN in this case,
    # propagating to the helium and total densities.
    @test AtmosphericModels._jr1971_helium_seasonal_correction(deg2rad(45), 0.0) == 0.0

    # The value must match the original formulation when the Sun declination is not zero.
    for δs in (-0.4, -0.1, 0.2), ϕ_gd in (-1.0, 0.3, 1.2)
        expected =
            0.65 / deg2rad(23.439291) *
            abs(δs) *
            (sin(π / 4 - ϕ_gd * δs / (2abs(δs)))^3 - 0.35355)

        @test AtmosphericModels._jr1971_helium_seasonal_correction(ϕ_gd, δs) == expected
    end
end

@testset "Errors" begin
    @test_throws ArgumentError AtmosphericModels.jr1971(now(), 0, 0, 89.9e3, 100, 100, 3)
    @test_throws ArgumentError AtmosphericModels.jr1971(now(), 0, 0, 3000.1e3, 100, 100, 3)
    @test_throws ArgumentError AtmosphericModels.jr1971(now(), 0, 0, NaN, 100, 100, 3)
end

@testset "Integer Inputs" begin
    # All-integer inputs must promote to a floating-point output instead of throwing an
    # InexactError when constructing the output structure.
    result   = AtmosphericModels.jr1971(2460000, 0, 0, 300_000, 100, 100, 3)
    expected = AtmosphericModels.jr1971(2460000.0, 0.0, 0.0, 300e3, 100.0, 100.0, 3.0)
    @test result isa AtmosphericModels.JR1971Output{Float64}
    @test result.total_density == expected.total_density
end

@testset "Float32 Output Type" begin
    # The output element type must be the promotion of the input types, and the call must
    # be type stable.
    result = AtmosphericModels.jr1971(
        2460000.25f0, 0.5f0, 0.5f0, 300.0f3, 100.0f0, 100.0f0, 3.0f0
    )
    @test result isa AtmosphericModels.JR1971Output{Float32}
end
