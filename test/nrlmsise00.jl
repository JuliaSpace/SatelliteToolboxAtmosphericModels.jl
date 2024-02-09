## Description #############################################################################
#
# Tests related to the NRLMSISE-00 Atmospheric Model.
#
## References ##############################################################################
#
# [1] https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/nrlmsise00_output.txt
# [2] https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
#
############################################################################################

# == Function: nrlmsise00 ==================================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# == Scenario 01 ===========================================================================
#
# Simulation using only the daily AP.
#
# Data obtained from the online version of NRLMSISE-00 [2]:
#
# Input parameters
#
# year= 1986, month= 6, day= 19, hour=21.50,
# Time_type = Universal
# Coordinate_type = Geographic
# latitude= -16.00, longitude= 312.00, height= 100.00
# Prof. parameters: start= 0.00 stop= 1000.00 step= 100.00
#
# Optional parametes: F10.7(daily) =121.; F10.7(3-month avg) =80.; ap(daily) = 7.
#
#    Selected parameters are:
# 1 Height, km
# 2 O, cm-3
# 3 N2, cm-3
# 4 O2, cm-3
# 5 Mass_density, g/cm-3
# 6 Temperature_neutral, K
# 7 Temperature_exospheric, K
# 8 He, cm-3
# 9 Ar, cm-3
# 10 H, cm-3
# 11 N, cm-3
# 12 Anomalous_Oxygen, cm-3
#
#       1          2          3          4         5      6     7          8          9         10         11         12
#     0.0  0.000E+00  1.918E+19  5.145E+18 1.180E-03  297.7  1027  1.287E+14  2.294E+17  0.000E+00  0.000E+00  0.000E+00
#   100.0  4.244E+11  9.498E+12  2.240E+12 5.783E-10  165.9  1027  1.061E+08  9.883E+10  2.209E+07  3.670E+05  0.000E+00
#   200.0  2.636E+09  2.248E+09  1.590E+08 1.838E-13  829.3   909  5.491E+06  1.829E+06  2.365E+05  3.125E+07  1.802E-09
#   300.0  3.321E+08  6.393E+07  2.881E+06 1.212E-14  900.6   909  3.173E+06  1.156E+04  1.717E+05  6.708E+06  4.938E+00
#   400.0  5.088E+07  2.414E+06  6.838E+04 1.510E-15  907.7   909  1.979E+06  1.074E+02  1.518E+05  1.274E+06  1.237E+03
#   500.0  8.340E+06  1.020E+05  1.839E+03 2.410E-16  908.5   909  1.259E+06  1.170E+00  1.355E+05  2.608E+05  4.047E+03
#   600.0  1.442E+06  4.729E+03  5.500E+01 4.543E-17  908.6   909  8.119E+05  1.455E-02  1.214E+05  5.617E+04  4.167E+03
#   700.0  2.622E+05  2.394E+02  1.817E+00 1.097E-17  908.6   909  5.301E+05  2.050E-04  1.092E+05  1.264E+04  3.173E+03
#   800.0  4.999E+04  1.317E+01  6.608E-02 3.887E-18  908.6   909  3.503E+05  3.254E-06  9.843E+04  2.964E+03  2.246E+03
#   900.0  9.980E+03  7.852E-01  2.633E-03 1.984E-18  908.6   909  2.342E+05  5.794E-08  8.900E+04  7.237E+02  1.571E+03
#  1000.0  2.082E+03  5.055E-02  1.146E-04 1.244E-18  908.6   909  1.582E+05  1.151E-09  8.069E+04  1.836E+02  1.103E+03
#
# == Scenario 02 ===========================================================================
#
# Simulation using the AP vector.
#
# Data obtained from the online version of NRLMSISE-00 [2]:
#
# Input parameters
# year= 2016, month= 6, day= 1, hour=11.00,
# Time_type = Universal
# Coordinate_type = Geographic
# latitude= -23.00, longitude= 315.00, height= 100.00
# Prof. parameters: start= 0.00 stop= 1000.00 step= 100.00
#
# Optional parametes: F10.7(daily) =not specified; F10.7(3-month avg) =not specified; ap(daily) = not specified
#
#    Selected parameters are:
# 1 Height, km
# 2 O, cm-3
# 3 N2, cm-3
# 4 O2, cm-3
# 5 Mass_density, g/cm-3
# 6 Temperature_neutral, K
# 7 Temperature_exospheric, K
# 8 He, cm-3
# 9 Ar, cm-3
# 10 H, cm-3
# 11 N, cm-3
# 12 Anomalous_Oxygen, cm-3
# 13 F10_7_daily
# 14 F10_7_3_month_avg
# 15 ap_daily
# 16 ap_00_03_hr_prior
# 17 ap_03_06_hr_prior
# 18 ap_06_09_hr_prior
# 19 ap_09_12_hr_prior
# 20 ap_12_33_hr_prior
# 21 ap_33_59_hr_prior
#
#       1          2          3          4         5      6     7          8          9         10         11         12      13     14     15      16     17     18      19     20     21
#     0.0  0.000E+00  1.919E+19  5.149E+18 1.181E-03  296.0  1027  1.288E+14  2.296E+17  0.000E+00  0.000E+00  0.000E+00   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   100.0  6.009E+11  9.392E+12  2.213E+12 5.764E-10  165.2  1027  1.339E+08  9.580E+10  2.654E+07  2.737E+05  0.000E+00   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   200.0  3.882E+09  1.978E+09  1.341E+08 2.028E-13  687.5   706  2.048E+07  9.719E+05  5.172E+05  1.632E+07  1.399E-09   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   300.0  3.140E+08  2.481E+07  9.340E+05 9.666E-15  705.5   706  1.082E+07  1.879E+03  3.860E+05  2.205E+06  3.876E+00   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   400.0  2.855E+07  3.737E+05  7.743E+03 8.221E-16  706.0   706  5.938E+06  4.688E+00  3.317E+05  2.633E+05  9.738E+02   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   500.0  2.787E+06  6.372E+03  7.382E+01 9.764E-17  706.0   706  3.319E+06  1.397E-02  2.868E+05  3.424E+04  3.186E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   600.0  2.910E+05  1.222E+02  8.047E-01 2.079E-17  706.0   706  1.887E+06  4.919E-05  2.490E+05  4.742E+03  3.281E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   700.0  3.240E+04  2.622E+00  9.973E-03 8.474E-18  706.0   706  1.090E+06  2.034E-07  2.171E+05  6.946E+02  2.498E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   800.0  3.835E+03  6.264E-02  1.398E-04 4.665E-18  706.0   706  6.393E+05  9.809E-10  1.900E+05  1.074E+02  1.768E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   900.0  4.816E+02  1.659E-03  2.204E-06 2.817E-18  706.0   706  3.806E+05  5.481E-12  1.669E+05  1.747E+01  1.236E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#  1000.0  6.400E+01  4.853E-05  3.891E-08 1.772E-18  706.0   706  2.298E+05  3.528E-14  1.471E+05  2.988E+00  8.675E+02   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#
############################################################################################

@testset "Providing All Space Indices" begin
    # == Scenario 01 =======================================================================

    # Test outputs.
    expected = [
           0.0  0.000E+00  1.918E+19  5.145E+18 1.180E-03  297.7  1027  1.287E+14 2.294E+17  0.000E+00  0.000E+00  0.000E+00
         100.0  4.244E+11  9.498E+12  2.240E+12 5.783E-10  165.9  1027  1.061E+08 9.883E+10  2.209E+07  3.670E+05  0.000E+00
         200.0  2.636E+09  2.248E+09  1.590E+08 1.838E-13  829.3   909  5.491E+06 1.829E+06  2.365E+05  3.125E+07  1.802E-09
         300.0  3.321E+08  6.393E+07  2.881E+06 1.212E-14  900.6   909  3.173E+06 1.156E+04  1.717E+05  6.708E+06  4.938E+00
         400.0  5.088E+07  2.414E+06  6.838E+04 1.510E-15  907.7   909  1.979E+06 1.074E+02  1.518E+05  1.274E+06  1.237E+03
         500.0  8.340E+06  1.020E+05  1.839E+03 2.410E-16  908.5   909  1.259E+06 1.170E+00  1.355E+05  2.608E+05  4.047E+03
         600.0  1.442E+06  4.729E+03  5.500E+01 4.543E-17  908.6   909  8.119E+05 1.455E-02  1.214E+05  5.617E+04  4.167E+03
         700.0  2.622E+05  2.394E+02  1.817E+00 1.097E-17  908.6   909  5.301E+05 2.050E-04  1.092E+05  1.264E+04  3.173E+03
         800.0  4.999E+04  1.317E+01  6.608E-02 3.887E-18  908.6   909  3.503E+05 3.254E-06  9.843E+04  2.964E+03  2.246E+03
         900.0  9.980E+03  7.852E-01  2.633E-03 1.984E-18  908.6   909  2.342E+05 5.794E-08  8.900E+04  7.237E+02  1.571E+03
        1000.0  2.082E+03  5.055E-02  1.146E-04 1.244E-18  908.6   909  1.582E+05 1.151E-09  8.069E+04  1.836E+02  1.103E+03
    ]

    # Constant input parameters.
    year    = 1986
    month   = 6
    day     = 19
    hour    = 21
    minute  = 30
    second  = 00
    ϕ_gd    = -16 * π / 180
    λ       = 312 * π / 180
    F10     = 121
    F10ₐ    = 80
    ap      = 7
    instant = date_to_jd(year, month, day, hour, minute, second) |> julian2datetime

    for i in 1:size(expected, 1)
        # Run the NRLMSISE-00 model wih the input parameters.
        out = AtmosphericModels.nrlmsise00(
            instant,
            expected[i, 1] * 1000,
            ϕ_gd,
            λ,
            F10ₐ,
            F10,
            ap;
            include_anomalous_oxygen = false
        )

        @test out.total_density           ≈ (expected[i,  5] * 1e3) rtol = 1e-3
        @test out.temperature             ≈ (expected[i,  6]      ) rtol = 1e-1 atol = 1e-9
        @test out.exospheric_temperature  ≈ (expected[i,  7]      ) rtol = 1e-1 atol = 1e-9
        @test out.O_number_density        ≈ (expected[i,  2] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N2_number_density       ≈ (expected[i,  3] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.O2_number_density       ≈ (expected[i,  4] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.He_number_density       ≈ (expected[i,  8] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.Ar_number_density       ≈ (expected[i,  9] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.H_number_density        ≈ (expected[i, 10] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N_number_density        ≈ (expected[i, 11] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.aO_number_density       ≈ (expected[i, 12] * 1e6) rtol = 1e-3 atol = 1e-9

        # Test the version that calls `gtd7d` instead of `gtd7`.
        out = AtmosphericModels.nrlmsise00(
            instant,
            expected[i, 1] * 1000,
            ϕ_gd,
            λ,
            F10ₐ,
            F10,
            ap;
            include_anomalous_oxygen = true
        )

        @test out.temperature             ≈ (expected[i,  6]      ) rtol = 1e-1 atol = 1e-9
        @test out.exospheric_temperature  ≈ (expected[i,  7]      ) rtol = 1e-1 atol = 1e-9
        @test out.O_number_density        ≈ (expected[i,  2] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N2_number_density       ≈ (expected[i,  3] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.O2_number_density       ≈ (expected[i,  4] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.He_number_density       ≈ (expected[i,  8] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.Ar_number_density       ≈ (expected[i,  9] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.H_number_density        ≈ (expected[i, 10] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N_number_density        ≈ (expected[i, 11] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.aO_number_density       ≈ (expected[i, 12] * 1e6) rtol = 1e-3 atol = 1e-9

        expected_total_density = expected[i, 5] + 1.66e-24 * 16 * expected[i, 12]
        @test out.total_density ≈ (expected_total_density * 1e3) rtol = 1e-3
    end

    # == Scenario 02 =======================================================================

    # Test outputs.
    expected = [
           0.0  0.000E+00  1.919E+19  5.149E+18 1.181E-03  296.0  1027  1.288E+14 2.296E+17  0.000E+00  0.000E+00  0.000E+00   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         100.0  6.009E+11  9.392E+12  2.213E+12 5.764E-10  165.2  1027  1.339E+08 9.580E+10  2.654E+07  2.737E+05  0.000E+00   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         200.0  3.882E+09  1.978E+09  1.341E+08 2.028E-13  687.5   706  2.048E+07 9.719E+05  5.172E+05  1.632E+07  1.399E-09   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         300.0  3.140E+08  2.481E+07  9.340E+05 9.666E-15  705.5   706  1.082E+07 1.879E+03  3.860E+05  2.205E+06  3.876E+00   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         400.0  2.855E+07  3.737E+05  7.743E+03 8.221E-16  706.0   706  5.938E+06 4.688E+00  3.317E+05  2.633E+05  9.738E+02   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         500.0  2.787E+06  6.372E+03  7.382E+01 9.764E-17  706.0   706  3.319E+06 1.397E-02  2.868E+05  3.424E+04  3.186E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         600.0  2.910E+05  1.222E+02  8.047E-01 2.079E-17  706.0   706  1.887E+06 4.919E-05  2.490E+05  4.742E+03  3.281E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         700.0  3.240E+04  2.622E+00  9.973E-03 8.474E-18  706.0   706  1.090E+06 2.034E-07  2.171E+05  6.946E+02  2.498E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         800.0  3.835E+03  6.264E-02  1.398E-04 4.665E-18  706.0   706  6.393E+05 9.809E-10  1.900E+05  1.074E+02  1.768E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
         900.0  4.816E+02  1.659E-03  2.204E-06 2.817E-18  706.0   706  3.806E+05 5.481E-12  1.669E+05  1.747E+01  1.236E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
        1000.0  6.400E+01  4.853E-05  3.891E-08 1.772E-18  706.0   706  2.298E+05 3.528E-14  1.471E+05  2.988E+00  8.675E+02   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6
    ]

    # Constant input parameters.
    year    = 2016
    month   = 6
    day     = 1
    hour    = 11
    minute  = 00
    second  = 00
    ϕ_gd    = -23 * π / 180
    λ       = -45 * π / 180
    instant = date_to_jd(year, month, day, hour, minute, second) |> julian2datetime

    # TODO: The tolerances here are much higher than in the previous test with the daily AP
    # only. This must be analyzed. However, it was observed that in the NRLMSISE-00 version,
    # the last component of the AP array is the "Average of eight 3 hour AP indices from 36
    # to 57 hours prior to current time." Whereas the label from the online version is
    # "ap_33_59_hr_prior", which seems to be different. On the other hand, the MSISE90
    # describe the last vector of AP array as "Average of eight 3 hr AP indicies [sic] from
    # 36 to 59 hrs prior to current time." Hence, it seems that the online version is using
    # the old algorithm related to the AP array. We must change the algorithm in our
    # implementation to the old one and see if the result is more similar to the online
    # version.
    #
    # The information can be obtained at:
    #
    #   https://ccmc.gsfc.nasa.gov/modelweb/
    #   https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/msise90/
    #

    for i in 1:size(expected, 1)
        F10  = expected[i, 13]
        F10ₐ = expected[i, 14]
        ap_a = expected[i, 15:21]

        # Run the NRLMSISE-00 model wih the input parameters.
        out = AtmosphericModels.nrlmsise00(
            instant,
            expected[i, 1] * 1000,
            ϕ_gd,
            λ,
            F10ₐ,
            F10,
            ap_a;
            include_anomalous_oxygen = false
        )

        @test out.total_density           ≈ (expected[i,  5] * 1e3) rtol = 5e-2
        @test out.temperature             ≈ (expected[i,  6]      ) rtol = 1e-1 atol = 1e-9
        @test out.exospheric_temperature  ≈ (expected[i,  7]      ) rtol = 1e-1 atol = 1e-9
        @test out.O_number_density        ≈ (expected[i,  2] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.N2_number_density       ≈ (expected[i,  3] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.O2_number_density       ≈ (expected[i,  4] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.He_number_density       ≈ (expected[i,  8] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.Ar_number_density       ≈ (expected[i,  9] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.H_number_density        ≈ (expected[i, 10] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.N_number_density        ≈ (expected[i, 11] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.aO_number_density       ≈ (expected[i, 12] * 1e6) rtol = 9e-2 atol = 1e-9

        # Test the version that calls `gtd7d` instead of `gtd7`.
        out = AtmosphericModels.nrlmsise00(
            instant,
            expected[i, 1] * 1000,
            ϕ_gd,
            λ,
            F10ₐ,
            F10,
            ap_a;
            include_anomalous_oxygen = true
        )

        @test out.temperature             ≈ (expected[i,  6]      ) rtol = 1e-1 atol = 1e-9
        @test out.exospheric_temperature  ≈ (expected[i,  7]      ) rtol = 1e-1 atol = 1e-9
        @test out.O_number_density        ≈ (expected[i,  2] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.N2_number_density       ≈ (expected[i,  3] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.O2_number_density       ≈ (expected[i,  4] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.He_number_density       ≈ (expected[i,  8] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.Ar_number_density       ≈ (expected[i,  9] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.H_number_density        ≈ (expected[i, 10] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.N_number_density        ≈ (expected[i, 11] * 1e6) rtol = 5e-2 atol = 1e-9
        @test out.aO_number_density       ≈ (expected[i, 12] * 1e6) rtol = 9e-2 atol = 1e-9

        expected_total_density = expected[i, 5] + 1.66e-24 * 16 * expected[i, 12]
        @test out.total_density ≈ (expected_total_density * 1e3) rtol = 5e-2
    end
end

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# == Scenario 01 ===========================================================================
#
# Data obtained for the instant 2023-01-01T10:00:00.000 using the online version of
# NRLMSISE-00 [2].
#
# Input parameters
# year= 2023, month= 1, day= 1, hour=10.00,
# Time_type = Universal
# Coordinate_type = Geographic
# latitude= -23.00, longitude= 315.00, height= 100.00
# Prof. parameters: start= 0.00 stop= 1000.00 step= 50.00
# Optional parametes: F10.7(daily) =not specified; F10.7(3-month avg) =not specified; ap(daily) = not specified
#
#    Selected parameters are:
# 1 Height, km
# 2 O, cm-3
# 3 N2, cm-3
# 4 O2, cm-3
# 5 Mass_density, g/cm-3
# 6 Temperature_neutral, K
# 7 Temperature_exospheric, K
# 8 He, cm-3
# 9 Ar, cm-3
# 10 H, cm-3
# 11 N, cm-3
# 12 Anomalous_Oxygen, cm-3
# 13 F10_7_daily
# 14 F10_7_3_month_avg
# 15 ap_daily
# 16 ap_00_03_hr_prior
# 17 ap_03_06_hr_prior
# 18 ap_06_09_hr_prior
# 19 ap_09_12_hr_prior
# 20 ap_12_33_hr_prior
# 21 ap_33_59_hr_prior
#
#       1          2          3          4         5      6     7          8          9         10         11         12      13     14     15      16     17     18      19     20     21
#     0.0  0.000E+00  1.909E+19  5.122E+18 1.175E-03  297.5  1027  1.281E+14  2.284E+17  0.000E+00  0.000E+00  0.000E+00  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#    50.0  0.000E+00  1.777E+16  4.766E+15 1.093E-06  269.0  1027  1.192E+11  2.125E+14  0.000E+00  0.000E+00  0.000E+00  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   100.0  3.854E+11  7.463E+12  1.511E+12 4.423E-10  188.3  1027  8.858E+07  7.379E+10  1.467E+07  2.513E+05  0.000E+00  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   150.0  1.703E+10  3.261E+10  1.970E+09 2.077E-12  688.5  1037  1.166E+07  6.061E+07  3.954E+05  1.246E+07  6.083E-21  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   200.0  4.318E+09  3.615E+09  1.381E+08 2.909E-13  911.5  1037  7.397E+06  2.827E+06  1.103E+05  2.694E+07  8.673E-09  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   250.0  1.584E+09  6.673E+08  1.955E+07 7.452E-14  991.0  1037  5.578E+06  2.589E+05  8.532E+04  1.422E+07  1.542E-02  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   300.0  6.521E+08  1.443E+08  3.382E+06 2.439E-14 1019.9  1037  4.421E+06  2.926E+04  7.816E+04  6.497E+06  2.402E+01  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   350.0  2.812E+08  3.338E+07  6.343E+05 9.150E-15 1030.5  1037  3.569E+06  3.624E+03  7.359E+04  3.058E+06  9.777E+02  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   400.0  1.243E+08  8.016E+06  1.243E+05 3.733E-15 1034.4  1037  2.906E+06  4.726E+02  6.977E+04  1.485E+06  6.035E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   450.0  5.580E+07  1.977E+06  2.509E+04 1.608E-15 1035.9  1037  2.377E+06  6.401E+01  6.631E+04  7.348E+05  1.404E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   500.0  2.540E+07  4.988E+05  5.200E+03 7.196E-16 1036.5  1037  1.952E+06  8.950E+00  6.311E+04  3.685E+05  1.975E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   550.0  1.170E+07  1.285E+05  1.104E+03 3.319E-16 1036.7  1037  1.608E+06  1.289E+00  6.012E+04  1.870E+05  2.142E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   600.0  5.450E+06  3.376E+04  2.396E+02 1.575E-16 1036.8  1037  1.329E+06  1.911E-01  5.731E+04  9.586E+04  2.034E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   650.0  2.567E+06  9.044E+03  5.317E+01 7.717E-17 1036.9  1037  1.101E+06  2.910E-02  5.468E+04  4.961E+04  1.805E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   700.0  1.222E+06  2.468E+03  1.205E+01 3.934E-17 1036.9  1037  9.143E+05  4.553E-03  5.220E+04  2.592E+04  1.548E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   750.0  5.882E+05  6.862E+02  2.791E+00 2.111E-17 1036.9  1037  7.615E+05  7.311E-04  4.987E+04  1.367E+04  1.307E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   800.0  2.859E+05  1.942E+02  6.595E-01 1.207E-17 1036.9  1037  6.358E+05  1.205E-04  4.767E+04  7.270E+03  1.096E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   850.0  1.404E+05  5.593E+01  1.590E-01 7.432E-18 1036.9  1037  5.323E+05  2.035E-05  4.560E+04  3.902E+03  9.162E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   900.0  6.962E+04  1.639E+01  3.910E-02 4.937E-18 1036.9  1037  4.467E+05  3.524E-06  4.364E+04  2.112E+03  7.661E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#   950.0  3.486E+04  4.884E+00  9.801E-03 3.517E-18 1036.9  1037  3.757E+05  6.250E-07  4.179E+04  1.153E+03  6.412E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#  1000.0  1.762E+04  1.480E+00  2.504E-03 2.653E-18 1036.9  1037  3.168E+05  1.135E-07  4.005E+04  6.346E+02  5.377E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
#
############################################################################################

@testset "Fetching Space Indices Automatically" begin
    SpaceIndices.init()

    expected = [
           0.0  0.000E+00  1.909E+19  5.122E+18 1.175E-03  297.5  1027  1.281E+14  2.284E+17  0.000E+00  0.000E+00  0.000E+00  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
          50.0  0.000E+00  1.777E+16  4.766E+15 1.093E-06  269.0  1027  1.192E+11  2.125E+14  0.000E+00  0.000E+00  0.000E+00  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         100.0  3.854E+11  7.463E+12  1.511E+12 4.423E-10  188.3  1027  8.858E+07  7.379E+10  1.467E+07  2.513E+05  0.000E+00  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         150.0  1.703E+10  3.261E+10  1.970E+09 2.077E-12  688.5  1037  1.166E+07  6.061E+07  3.954E+05  1.246E+07  6.083E-21  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         200.0  4.318E+09  3.615E+09  1.381E+08 2.909E-13  911.5  1037  7.397E+06  2.827E+06  1.103E+05  2.694E+07  8.673E-09  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         250.0  1.584E+09  6.673E+08  1.955E+07 7.452E-14  991.0  1037  5.578E+06  2.589E+05  8.532E+04  1.422E+07  1.542E-02  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         300.0  6.521E+08  1.443E+08  3.382E+06 2.439E-14 1019.9  1037  4.421E+06  2.926E+04  7.816E+04  6.497E+06  2.402E+01  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         350.0  2.812E+08  3.338E+07  6.343E+05 9.150E-15 1030.5  1037  3.569E+06  3.624E+03  7.359E+04  3.058E+06  9.777E+02  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         400.0  1.243E+08  8.016E+06  1.243E+05 3.733E-15 1034.4  1037  2.906E+06  4.726E+02  6.977E+04  1.485E+06  6.035E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         450.0  5.580E+07  1.977E+06  2.509E+04 1.608E-15 1035.9  1037  2.377E+06  6.401E+01  6.631E+04  7.348E+05  1.404E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         500.0  2.540E+07  4.988E+05  5.200E+03 7.196E-16 1036.5  1037  1.952E+06  8.950E+00  6.311E+04  3.685E+05  1.975E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         550.0  1.170E+07  1.285E+05  1.104E+03 3.319E-16 1036.7  1037  1.608E+06  1.289E+00  6.012E+04  1.870E+05  2.142E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         600.0  5.450E+06  3.376E+04  2.396E+02 1.575E-16 1036.8  1037  1.329E+06  1.911E-01  5.731E+04  9.586E+04  2.034E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         650.0  2.567E+06  9.044E+03  5.317E+01 7.717E-17 1036.9  1037  1.101E+06  2.910E-02  5.468E+04  4.961E+04  1.805E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         700.0  1.222E+06  2.468E+03  1.205E+01 3.934E-17 1036.9  1037  9.143E+05  4.553E-03  5.220E+04  2.592E+04  1.548E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         750.0  5.882E+05  6.862E+02  2.791E+00 2.111E-17 1036.9  1037  7.615E+05  7.311E-04  4.987E+04  1.367E+04  1.307E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         800.0  2.859E+05  1.942E+02  6.595E-01 1.207E-17 1036.9  1037  6.358E+05  1.205E-04  4.767E+04  7.270E+03  1.096E+04  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         850.0  1.404E+05  5.593E+01  1.590E-01 7.432E-18 1036.9  1037  5.323E+05  2.035E-05  4.560E+04  3.902E+03  9.162E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         900.0  6.962E+04  1.639E+01  3.910E-02 4.937E-18 1036.9  1037  4.467E+05  3.524E-06  4.364E+04  2.112E+03  7.661E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
         950.0  3.486E+04  4.884E+00  9.801E-03 3.517E-18 1036.9  1037  3.757E+05  6.250E-07  4.179E+04  1.153E+03  6.412E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
        1000.0  1.762E+04  1.480E+00  2.504E-03 2.653E-18 1036.9  1037  3.168E+05  1.135E-07  4.005E+04  6.346E+02  5.377E+03  159.5  153.6   13.8   12.0    7.0   15.0    9.0   15.9   28.6
    ]

    for i in 1:size(expected, 1)
        # Run the NRLMSISE-00 model wih the input parameters.
        out = AtmosphericModels.nrlmsise00(
            DateTime("2023-01-01T10:00:00"),
            expected[i, 1] * 1000,
            -23 |> deg2rad,
            -45 |> deg2rad;
            include_anomalous_oxygen = false
        )

        @test out.total_density           ≈ (expected[i,  5] * 1e3) rtol = 3e-3
        @test out.temperature             ≈ (expected[i,  6]      ) rtol = 1e-1 atol = 1e-9
        @test out.exospheric_temperature  ≈ (expected[i,  7]      ) rtol = 1e-1 atol = 1e-9
        @test out.O_number_density        ≈ (expected[i,  2] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N2_number_density       ≈ (expected[i,  3] * 1e6) rtol = 3e-3 atol = 1e-9
        @test out.O2_number_density       ≈ (expected[i,  4] * 1e6) rtol = 3e-3 atol = 1e-9
        @test out.He_number_density       ≈ (expected[i,  8] * 1e6) rtol = 2e-3 atol = 1e-9
        @test out.Ar_number_density       ≈ (expected[i,  9] * 1e6) rtol = 3e-3 atol = 1e-9
        @test out.H_number_density        ≈ (expected[i, 10] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N_number_density        ≈ (expected[i, 11] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.aO_number_density       ≈ (expected[i, 12] * 1e6) rtol = 1e-3 atol = 1e-9
    end

    # For altitudes lower than 80 km, we use default space indices.
    expected = AtmosphericModels.nrlmsise00(
        DateTime("2023-01-01T10:00:00"),
        79e3,
        -23 |> deg2rad,
        -45 |> deg2rad,
        150,
        150,
        4
    )

    result = AtmosphericModels.nrlmsise00(
        DateTime("2023-01-01T10:00:00"),
        79e3,
        -23 |> deg2rad,
        -45 |> deg2rad
    )

    @test result.total_density           ≈ expected.total_density
    @test result.temperature             ≈ expected.temperature
    @test result.exospheric_temperature  ≈ expected.exospheric_temperature
    @test result.O_number_density        ≈ expected.O_number_density
    @test result.N2_number_density       ≈ expected.N2_number_density
    @test result.O2_number_density       ≈ expected.O2_number_density
    @test result.He_number_density       ≈ expected.He_number_density
    @test result.Ar_number_density       ≈ expected.Ar_number_density
    @test result.H_number_density        ≈ expected.H_number_density
    @test result.N_number_density        ≈ expected.N_number_density
    @test result.aO_number_density       ≈ expected.aO_number_density
end

@testset "Pre-allocating the Legendre Matrix" begin
    # We will use the scenario 1 of the first test.

    # Test outputs.
    expected = [
           0.0  0.000E+00  1.918E+19  5.145E+18 1.180E-03  297.7  1027  1.287E+14 2.294E+17  0.000E+00  0.000E+00  0.000E+00
         100.0  4.244E+11  9.498E+12  2.240E+12 5.783E-10  165.9  1027  1.061E+08 9.883E+10  2.209E+07  3.670E+05  0.000E+00
         200.0  2.636E+09  2.248E+09  1.590E+08 1.838E-13  829.3   909  5.491E+06 1.829E+06  2.365E+05  3.125E+07  1.802E-09
         300.0  3.321E+08  6.393E+07  2.881E+06 1.212E-14  900.6   909  3.173E+06 1.156E+04  1.717E+05  6.708E+06  4.938E+00
         400.0  5.088E+07  2.414E+06  6.838E+04 1.510E-15  907.7   909  1.979E+06 1.074E+02  1.518E+05  1.274E+06  1.237E+03
         500.0  8.340E+06  1.020E+05  1.839E+03 2.410E-16  908.5   909  1.259E+06 1.170E+00  1.355E+05  2.608E+05  4.047E+03
         600.0  1.442E+06  4.729E+03  5.500E+01 4.543E-17  908.6   909  8.119E+05 1.455E-02  1.214E+05  5.617E+04  4.167E+03
         700.0  2.622E+05  2.394E+02  1.817E+00 1.097E-17  908.6   909  5.301E+05 2.050E-04  1.092E+05  1.264E+04  3.173E+03
         800.0  4.999E+04  1.317E+01  6.608E-02 3.887E-18  908.6   909  3.503E+05 3.254E-06  9.843E+04  2.964E+03  2.246E+03
         900.0  9.980E+03  7.852E-01  2.633E-03 1.984E-18  908.6   909  2.342E+05 5.794E-08  8.900E+04  7.237E+02  1.571E+03
        1000.0  2.082E+03  5.055E-02  1.146E-04 1.244E-18  908.6   909  1.582E+05 1.151E-09  8.069E+04  1.836E+02  1.103E+03
    ]

    # Constant input parameters.
    year    = 1986
    month   = 6
    day     = 19
    hour    = 21
    minute  = 30
    second  = 00
    ϕ_gd    = -16 * π / 180
    λ       = 312 * π / 180
    F10     = 121
    F10ₐ    = 80
    ap      = 7
    instant = date_to_jd(year, month, day, hour, minute, second) |> julian2datetime
    P       = zeros(9, 9)

    for i in 1:size(expected, 1)
        # Initialize `P` to check if it is used inside `nrlmsise00`.
        P .= 0.0

        # Run the NRLMSISE-00 model wih the input parameters.
        out = AtmosphericModels.nrlmsise00(
            instant,
            expected[i, 1] * 1000,
            ϕ_gd,
            λ,
            F10ₐ,
            F10,
            ap;
            include_anomalous_oxygen = false,
            P = P
        )

        @test out.total_density           ≈ (expected[i,  5] * 1e3) rtol = 1e-3
        @test out.temperature             ≈ (expected[i,  6]      ) rtol = 1e-1 atol = 1e-9
        @test out.exospheric_temperature  ≈ (expected[i,  7]      ) rtol = 1e-1 atol = 1e-9
        @test out.O_number_density        ≈ (expected[i,  2] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N2_number_density       ≈ (expected[i,  3] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.O2_number_density       ≈ (expected[i,  4] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.He_number_density       ≈ (expected[i,  8] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.Ar_number_density       ≈ (expected[i,  9] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.H_number_density        ≈ (expected[i, 10] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N_number_density        ≈ (expected[i, 11] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.aO_number_density       ≈ (expected[i, 12] * 1e6) rtol = 1e-3 atol = 1e-9

        # Check that P was used.
        @test abs(sum(P)) > 0

        # Initialize `P` to check if it is used inside `nrlmsise00`.
        P .= 0.0

        # Test the version that calls `gtd7d` instead of `gtd7`.
        out = AtmosphericModels.nrlmsise00(
            instant,
            expected[i, 1] * 1000,
            ϕ_gd,
            λ,
            F10ₐ,
            F10,
            ap;
            include_anomalous_oxygen = true,
            P = P
        )

        @test out.temperature             ≈ (expected[i,  6]      ) rtol = 1e-1 atol = 1e-9
        @test out.exospheric_temperature  ≈ (expected[i,  7]      ) rtol = 1e-1 atol = 1e-9
        @test out.O_number_density        ≈ (expected[i,  2] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N2_number_density       ≈ (expected[i,  3] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.O2_number_density       ≈ (expected[i,  4] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.He_number_density       ≈ (expected[i,  8] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.Ar_number_density       ≈ (expected[i,  9] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.H_number_density        ≈ (expected[i, 10] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.N_number_density        ≈ (expected[i, 11] * 1e6) rtol = 1e-3 atol = 1e-9
        @test out.aO_number_density       ≈ (expected[i, 12] * 1e6) rtol = 1e-3 atol = 1e-9

        expected_total_density = expected[i, 5] + 1.66e-24 * 16 * expected[i, 12]
        @test out.total_density ≈ (expected_total_density * 1e3) rtol = 1e-3

        # Check that P was used.
        @test abs(sum(P)) > 0

    end
end

@testset "Errors" begin
    # == Wrong Size in Matrix `P` ==========================================================

    P = zeros(9, 8)
    @test_throws ArgumentError AtmosphericModels.nrlmsise00(
        now() |> datetime2julian,
        100e3,
        0,
        0,
        100,
        100,
        20;
        P = P
    )

    P = zeros(8, 9)
    @test_throws ArgumentError AtmosphericModels.nrlmsise00(
        now() |> datetime2julian,
        100e3,
        0,
        0,
        100,
        100,
        20;
        P = P
    )

    P = zeros(8, 8)
    @test_throws ArgumentError AtmosphericModels.nrlmsise00(
        now() |> datetime2julian,
        100e3,
        0,
        0,
        100,
        100,
        20;
        P = P
    )
end

@testset "Show" begin
    result = AtmosphericModels.nrlmsise00(
        DateTime("2023-01-01T10:00:00"),
        500e3,
        0,
        0,
        100,
        100,
        3
    )

    expected = "NRLMSISE-00 output (ρ = 2.5893e-13 kg / m³)"
    str = sprint(show, result)

    expected = """
        NRLMSISE-00 Atmospheric Model Result:
                  Total density :     2.5893e-13  kg / m³
                    Temperature :         875.50  K
               Exospheric Temp. :         875.56  K
              N  number density :    2.02467e+11  1 / m³
              N₂ number density :    7.80367e+10  1 / m³
              O  number density :      8.908e+12  1 / m³
          Anomalous O num. den. :    3.84174e+09  1 / m³
              O₂ number density :    1.05499e+09  1 / m³
              Ar number density :         619170  1 / m³
              He number density :    2.04514e+12  1 / m³
              H  number density :    1.58356e+11  1 / m³"""

    str = sprint(show, MIME("text/plain"), result)
    @test str == expected
end
