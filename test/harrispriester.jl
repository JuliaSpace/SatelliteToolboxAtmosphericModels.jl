## Description #############################################################################
#
# Tests related to the Harris-Priester models.
#
############################################################################################

# == Harris-Priester Standard Model ========================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# == Scenario 01 ===========================================================================
#
#   Values obtained from Orekit (Java) via orekit-jpype using the following inputs:
#
#       Date:       2003-03-21 01:00:00 UTC
#       Latitude:   0°
#       Longitude:  0°
#       n:          4 (default)
#
#   Result:
#
#       | Altitude [km] | Density [kg/m³]      |
#       |---------------|----------------------|
#       |           150 | 2.1220430185270577e-09|
#       |           200 | 2.5572798517082807e-10|
#       |           300 | 1.7088409428192282e-11|
#       |           400 | 2.2514252272834970e-12|
#       |           500 | 3.9236341695759752e-13|
#       |           600 | 8.0958249931790292e-14|
#       |           700 | 2.0521620211337461e-14|
#       |           800 | 7.1030082271810029e-15|
#       |           900 | 2.6748844701934813e-15|
#       |          1000 | 1.1578404734799309e-15|
#
# == Scenario 02a ==========================================================================
#
#   Same as Scenario 01 but with n = 2:
#
#       | Altitude [km] | Density [kg/m³]      |
#       |---------------|----------------------|
#       |           150 | 2.1240001807459317e-09|
#       |           200 | 2.5700119285084817e-10|
#       |           300 | 1.7471003074849554e-11|
#       |           400 | 2.3617628779668907e-12|
#       |           500 | 4.2709568067834383e-13|
#       |           600 | 9.2707536671545896e-14|
#       |           700 | 2.4689954842437929e-14|
#       |           800 | 8.6502396625984707e-15|
#       |           900 | 3.3828558672366406e-15|
#       |          1000 | 1.5145490714359715e-15|
#
# == Scenario 02b ==========================================================================
#
#   Same as Scenario 01 but with n = 6:
#
#       | Altitude [km] | Density [kg/m³]      |
#       |---------------|----------------------|
#       |           150 | 2.1220009252132209e-09|
#       |           200 | 2.5570060188602002e-10|
#       |           300 | 1.7080180864261531e-11|
#       |           400 | 2.2490521601388917e-12|
#       |           500 | 3.9161641905268500e-13|
#       |           600 | 8.0705554263883928e-14|
#       |           700 | 2.0431970505189843e-14|
#       |           800 | 7.0697314258194706e-15|
#       |           900 | 2.6596579029135220e-15|
#       |          1000 | 1.1501686275708982e-15|
#
# == Scenario 03 ===========================================================================
#
#   Values obtained from Orekit for latitude variation at 500 km:
#
#       Date:       2003-03-21 01:00:00 UTC
#       Longitude:  0°
#       Altitude:   500 km
#       n:          4
#
#       | Latitude [°] | Density [kg/m³]      |
#       |--------------|----------------------|
#       |            0 | 3.9236341695759752e-13|
#       |           15 | 3.9393220819088735e-13|
#       |           30 | 4.0350599682163622e-13|
#       |           45 | 4.3415897821552519e-13|
#       |           60 | 5.0285484446420926e-13|
#       |           75 | 6.2424630627311681e-13|
#       |          -30 | 4.0350008559554109e-13|
#       |          -60 | 5.0282346210873605e-13|
#
# == Scenario 04 ===========================================================================
#
#   Values obtained from Orekit for longitude (diurnal bulge) variation:
#
#       Date:       2003-03-21 01:00:00 UTC
#       Latitude:   0°
#       Altitude:   500 km
#       n:          4
#
#       | Longitude [°] | Density [kg/m³]      |
#       |---------------|----------------------|
#       |             0 | 3.9236341695759752e-13|
#       |            30 | 3.9188232351247705e-13|
#       |            60 | 4.2173248147961520e-13|
#       |            90 | 5.9951186530463867e-13|
#       |           120 | 1.0130114879877873e-12|
#       |           150 | 1.5612650147206855e-12|
#       |           180 | 1.9717720556009100e-12|
#       |           210 | 1.9991107252102569e-12|
#       |           240 | 1.6261250939445796e-12|
#       |           270 | 1.0783526606567886e-12|
#       |           300 | 6.3799546180942567e-13|
#       |           330 | 4.3287740837505535e-13|
#
# == Scenario 05 ===========================================================================
#
#   Values obtained from Orekit for a different epoch (summer solstice):
#
#       Date:       2023-06-21 12:00:00 UTC
#       Latitude:   0°
#       Longitude:  0°
#       n:          4
#
#       | Altitude [km] | Density [kg/m³]      |
#       |---------------|----------------------|
#       |           150 | 2.1965736409949631e-09|
#       |           200 | 3.0421296000227622e-10|
#       |           300 | 3.1657944013913403e-11|
#       |           400 | 6.4531892444965479e-12|
#       |           500 | 1.7150014741784003e-12|
#       |           600 | 5.2838240610384661e-13|
#       |           700 | 1.7925581797777100e-13|
#       |           800 | 6.6023071608745029e-14|
#       |           900 | 2.9635117200251611e-14|
#       |          1000 | 0.0000000000000000e+00|
#
############################################################################################

@testset "Providing All Parameters" begin
    # Common inputs.
    jd      = date_to_jd(2003, 3, 21, 1, 0, 0)
    instant = julian2datetime(jd)
    h       = [150, 200, 300, 400, 500, 600, 700, 800, 900, 1000] * 1000

    # == Scenario 01: Equator, lon=0, n=4 (default) =======================================

    ϕ_gd = 0.0
    λ    = 0.0

    # Results from Orekit [kg/m³].
    results_n4 = [
        2.1220430185270577e-09
        2.5572798517082807e-10
        1.7088409428192282e-11
        2.2514252272834970e-12
        3.9236341695759752e-13
        8.0958249931790292e-14
        2.0521620211337461e-14
        7.1030082271810029e-15
        2.6748844701934813e-15
        1.1578404734799309e-15
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester(instant, ϕ_gd, λ, h[i])
        @test ret ≈ results_n4[i] rtol=5e-4

        ret_jd = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h[i])
        @test ret_jd ≈ results_n4[i] rtol=5e-4
    end

    # == Scenario 02a: Equator, lon=0, n=2 ================================================

    results_n2 = [
        2.1240001807459317e-09
        2.5700119285084817e-10
        1.7471003074849554e-11
        2.3617628779668907e-12
        4.2709568067834383e-13
        9.2707536671545896e-14
        2.4689954842437929e-14
        8.6502396625984707e-15
        3.3828558672366406e-15
        1.5145490714359715e-15
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h[i]; n=2)
        @test ret ≈ results_n2[i] rtol=5e-4
    end

    # == Scenario 02b: Equator, lon=0, n=6 ================================================

    results_n6 = [
        2.1220009252132209e-09
        2.5570060188602002e-10
        1.7080180864261531e-11
        2.2490521601388917e-12
        3.9161641905268500e-13
        8.0705554263883928e-14
        2.0431970505189843e-14
        7.0697314258194706e-15
        2.6596579029135220e-15
        1.1501686275708982e-15
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h[i]; n=6)
        @test ret ≈ results_n6[i] rtol=5e-4
    end
end

@testset "Latitude Variation" begin
    jd   = date_to_jd(2003, 3, 21, 1, 0, 0)
    λ    = 0.0
    h    = 500e3

    # Results from Orekit [kg/m³].
    # Wider tolerance needed because latitude affects the Sun position angle calculation and
    # Orekit uses a different Sun ephemeris than SatelliteToolboxCelestialBodies.
    lats_deg = [0, 15, 30, 45, 60, 75, -30, -60]
    results  = [
        3.9236341695759752e-13
        3.9393220819088735e-13
        4.0350599682163622e-13
        4.3415897821552519e-13
        5.0285484446420926e-13
        6.2424630627311681e-13
        4.0350008559554109e-13
        5.0282346210873605e-13
    ]

    for i in eachindex(lats_deg)
        ϕ_gd = deg2rad(lats_deg[i])
        ret = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h)
        @test ret ≈ results[i] rtol=3e-3
    end
end

@testset "Diurnal Bulge (Longitude Variation)" begin
    jd   = date_to_jd(2003, 3, 21, 1, 0, 0)
    ϕ_gd = 0.0
    h    = 500e3

    # Results from Orekit [kg/m³].
    lons_deg = collect(0:30:330)
    results  = [
        3.9236341695759752e-13
        3.9188232351247705e-13
        4.2173248147961520e-13
        5.9951186530463867e-13
        1.0130114879877873e-12
        1.5612650147206855e-12
        1.9717720556009100e-12
        1.9991107252102569e-12
        1.6261250939445796e-12
        1.0783526606567886e-12
        6.3799546180942567e-13
        4.3287740837505535e-13
    ]

    for i in eachindex(lons_deg)
        λ = deg2rad(lons_deg[i])
        ret = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h)
        @test ret ≈ results[i] rtol=5e-4
    end
end

@testset "Different Epoch (Summer Solstice)" begin
    jd   = date_to_jd(2023, 6, 21, 12, 0, 0)
    ϕ_gd = 0.0
    λ    = 0.0
    h    = [150, 200, 300, 400, 500, 600, 700, 800, 900] * 1000

    # Results from Orekit [kg/m³].
    # NOTE: Orekit returns 0.0 at 1000 km for this epoch; we exclude that since the
    # density table max is 1000 km and our implementation also returns 0 above max.
    results = [
        2.1965736409949631e-09
        3.0421296000227622e-10
        3.1657944013913403e-11
        6.4531892444965479e-12
        1.7150014741784003e-12
        5.2838240610384661e-13
        1.7925581797777100e-13
        6.6023071608745029e-14
        2.9635117200251611e-14
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h[i])
        @test ret ≈ results[i] rtol=5e-4
    end
end

@testset "Boundary Conditions" begin
    jd   = date_to_jd(2003, 3, 21, 1, 0, 0)
    ϕ_gd = 0.0
    λ    = 0.0

    # Above maximum altitude returns 0.
    @test AtmosphericModels.harrispriester(jd, ϕ_gd, λ, 1500e3) == 0.0

    # Below minimum altitude throws.
    @test_throws ArgumentError AtmosphericModels.harrispriester(jd, ϕ_gd, λ, 50e3)
end

@testset "Custom Density Table" begin
    jd   = date_to_jd(2003, 3, 21, 1, 0, 0)
    ϕ_gd = 0.0
    λ    = 0.0

    user_tab = [
        100000.   4.974e+02  4.974e+02;
        110000.   7.800e+01  7.800e+01;
        120000.   2.490e+01  2.400e+01;
        130000.   8.377e+00  8.710e+00;
        140000.   3.899e+00  4.059e+00;
        150000.   2.122e+00  2.215e+00;
        160000.   1.263e+00  1.344e+00;
        170000.   8.008e-01  8.758e-01;
        180000.   5.283e-01  6.010e-01;
        190000.   3.618e-01  4.297e-01;
        200000.   2.557e-01  3.162e-01;
        210000.   1.839e-01  2.396e-01;
        220000.   1.341e-01  1.853e-01;
        230000.   9.949e-02  1.455e-01;
        240000.   7.488e-02  1.157e-01;
        250000.   5.709e-02  9.308e-02;
        260000.   4.403e-02  7.555e-02;
        270000.   3.430e-02  6.182e-02;
        280000.   2.697e-02  5.095e-02;
        290000.   2.139e-02  4.226e-02;
        300000.   1.708e-02  3.526e-02;
        320000.   1.099e-02  2.511e-02;
        340000.   7.214e-03  1.819e-02;
        360000.   4.824e-03  1.337e-02;
        380000.   3.274e-03  9.955e-03;
        400000.   2.249e-03  7.492e-03;
        420000.   1.558e-03  5.684e-03;
        440000.   1.091e-03  4.355e-03;
        460000.   7.701e-04  3.362e-03;
        480000.   5.474e-04  2.612e-03;
        500000.   3.916e-04  2.042e-03;
        520000.   2.819e-04  1.605e-03;
        540000.   2.042e-04  1.267e-03;
        560000.   1.488e-04  1.005e-03;
        580000.   1.092e-04  7.997e-04;
        600000.   8.070e-05  6.390e-04;
        620000.   6.012e-05  5.123e-04;
        640000.   4.519e-05  4.121e-04;
        660000.   3.430e-05  3.325e-04;
        680000.   2.632e-05  2.691e-04;
        700000.   2.043e-05  2.185e-04;
        720000.   1.607e-05  1.779e-04;
        740000.   1.281e-05  1.452e-04;
        760000.   1.036e-05  1.190e-04;
        780000.   8.496e-06  9.776e-05;
        800000.   7.069e-06  8.059e-05;
        850000.   4.800e-06  5.500e-05;
        900000.   3.300e-06  3.700e-05;
        950000.   2.450e-06  2.400e-05;
        1000000.  1.900e-06  1.700e-05;
        1100000.  1.180e-06  8.700e-06;
        1200000.  7.500e-07  4.800e-06;
        1300000.  5.300e-07  3.200e-06;
        1400000.  4.100e-07  2.000e-06;
        1500000.  2.900e-07  1.350e-06;
        1600000.  2.000e-07  9.500e-07;
        1700000.  1.600e-07  7.700e-07;
        1800000.  1.200e-07  6.300e-07;
        1900000.  9.600e-08  5.200e-07;
        2000000.  7.300e-08  4.400e-07
    ]

    # With the user table, 1500 km is within range (max is 2000 km) and should return > 0.
    rho = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, 1500e3; alt_ρ=user_tab)
    @test rho > 0

    # Above max altitude (2000 km) should return 0.
    rho_above = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, 2500e3; alt_ρ=user_tab)
    @test rho_above == 0.0
end

# == Harris-Priester Modified Model (Hatten & Russell 2017) ================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# All reference values obtained by compiling and running the original FORTRAN
# implementation (harris_priester_mod_dist.f90) by Noble Hatten and Ryan P. Russell.
#
# The FORTRAN code takes ECI position/velocity; the Julia implementation takes geodetic
# coordinates. The driver program converts geodetic to ECI for consistency. Small
# discrepancies arise from:
#   1. Different Sun position algorithms (USNO vs SatelliteToolboxCelestialBodies)
#   2. Different GMST computation methods
#
# The FORTRAN code computes n from orbital inclination:
#   n = 2.001 + 4*sin²(inclination)
# For equatorial orbits (i=0°): n = 2.001
#
# == Scenario 01 ===========================================================================
#
#   FORTRAN reference for equator (lat=0°, lon=0°) with varying altitude and F10.7:
#
#   JD = 2452720.541666667 (2003-03-21 01:00:00 UTC), equatorial orbit (n=2.001)
#
#       | Alt [km] | F10.7=75 [kg/m³]     | F10.7=150 [kg/m³]    | F10.7=250 [kg/m³]    |
#       |----------|----------------------|----------------------|----------------------|
#       |      120 | 1.648445e-08         | 1.700478e-08         | 1.766926e-08         |
#       |      200 | 1.790835e-10         | 2.631724e-10         | 3.702013e-10         |
#       |      300 | 7.981701e-12         | 2.177584e-11         | 4.383177e-11         |
#       |      400 | 6.527816e-13         | 2.994375e-12         | 7.981621e-12         |
#       |      500 | 7.629473e-14         | 5.260829e-13         | 1.753297e-12         |
#       |      600 | 1.415497e-14         | 1.128321e-13         | 4.423989e-13         |
#       |      700 | 4.546487e-15         | 3.047021e-14         | 1.286127e-13         |
#       |      800 | 2.248141e-15         | 1.057779e-14         | 4.407029e-14         |
#       |      900 | 1.385529e-15         | 4.680785e-15         | 1.812867e-14         |
#       |     1000 | 9.260778e-16         | 2.537412e-15         | 8.531580e-15         |
#
############################################################################################

@testset "Providing F10.7" begin
    jd      = date_to_jd(2003, 3, 21, 1, 0, 0)
    instant = julian2datetime(jd)
    ϕ_gd    = 0.0
    λ       = 0.0
    h       = [120, 200, 300, 400, 500, 600, 700, 800, 900, 1000] * 1000

    # n value for equatorial orbit (FORTRAN: npow = 2.001 + 4*sin²(0°) = 2.001)
    n_eq = 2.001

    # == Low Solar Activity: F10.7 = 75 ===================================================

    results_f75 = [
        1.6484454506680808e-08
        1.7908348285099853e-10
        7.9817012236982133e-12
        6.5278163368528728e-13
        7.6294733752142897e-14
        1.4154969529780225e-14
        4.5464870830758152e-15
        2.2481413326961890e-15
        1.3855288566966760e-15
        9.2607782765970119e-16
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester_modified(instant, ϕ_gd, λ, h[i], 75.0; n=n_eq)
        @test ret ≈ results_f75[i] rtol=2e-2

        ret_jd = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, h[i], 75.0; n=n_eq)
        @test ret_jd ≈ results_f75[i] rtol=2e-2
    end

    # == Medium Solar Activity: F10.7 = 150 ===============================================

    results_f150 = [
        1.7004779502658466e-08
        2.6317237403622035e-10
        2.1775835950869826e-11
        2.9943745351654492e-12
        5.2608287456757322e-13
        1.1283206282989865e-13
        3.0470213172162873e-14
        1.0577793770567949e-14
        4.6807852698593273e-15
        2.5374123357717461e-15
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester_modified(instant, ϕ_gd, λ, h[i], 150.0; n=n_eq)
        @test ret ≈ results_f150[i] rtol=2e-2
    end

    # == High Solar Activity: F10.7 = 250 =================================================

    results_f250 = [
        1.7669257863183996e-08
        3.7020131535920146e-10
        4.3831767874961432e-11
        7.9816213217197887e-12
        1.7532974268294977e-12
        4.4239891581501132e-13
        1.2861269369088979e-13
        4.4070285359495685e-14
        1.8128672970361650e-14
        8.5315804688572281e-15
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester_modified(instant, ϕ_gd, λ, h[i], 250.0; n=n_eq)
        @test ret ≈ results_f250[i] rtol=2e-2
    end
end

@testset "Latitude Variation (Modified)" begin
    # NOTE: The FORTRAN reference code computes n from orbital inclination internally
    # (npow = 2.001 + 4*sin²(i)), coupling n to the satellite position/velocity. The Julia
    # implementation takes n as an explicit parameter, so non-equatorial FORTRAN reference
    # values cannot be directly compared at a fixed n. We test physical behavior instead.
    jd   = date_to_jd(2003, 3, 21, 1, 0, 0)
    h    = 500e3
    F10ₐ = 150.0

    # All densities must be positive at various latitudes.
    for lat_deg in [-60, -30, 0, 30, 60]
        ϕ_gd = deg2rad(lat_deg)
        ret = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, 0.0, h, F10ₐ)
        @test ret > 0
    end

    # Symmetric latitudes should give similar (not necessarily identical) densities due
    # to the Sun's declination at vernal equinox (~0°).
    rho_n30 = AtmosphericModels.harrispriester_modified(jd, deg2rad(-30.0), 0.0, h, F10ₐ)
    rho_p30 = AtmosphericModels.harrispriester_modified(jd, deg2rad(30.0),  0.0, h, F10ₐ)
    @test rho_n30 ≈ rho_p30 rtol=0.1
end

@testset "Diurnal Bulge (Modified)" begin
    jd    = date_to_jd(2003, 3, 21, 1, 0, 0)
    h     = 500e3
    F10ₐ  = 150.0
    n_eq  = 2.001

    # FORTRAN reference values for longitude variation.
    lons_deg = collect(0:30:330)
    results  = [
        5.2608287456757322e-13
        5.1865112512137276e-13
        6.3206225678180850e-13
        8.3606594213847372e-13
        1.0760098841637155e-12
        1.2875733594800382e-12
        1.4140326704576779e-12
        1.4214794839272383e-12
        1.3079169964396878e-12
        1.1037955771451819e-12
        8.6384427868740631e-13
        6.5238797658846279e-13
    ]

    for i in eachindex(lons_deg)
        λ = deg2rad(lons_deg[i])
        ret = AtmosphericModels.harrispriester_modified(jd, 0.0, λ, h, F10ₐ; n=n_eq)
        @test ret ≈ results[i] rtol=2e-2
    end
end

@testset "Different Epoch (Modified)" begin
    jd   = date_to_jd(2023, 6, 21, 12, 0, 0)
    ϕ_gd = 0.0
    λ    = 0.0
    h    = [120, 200, 300, 400, 500, 600, 700, 800, 900, 1000] * 1000
    F10ₐ = 150.0
    n_eq = 2.001

    # FORTRAN reference values.
    results = [
        2.3451895104567802e-08
        3.1955247601053275e-10
        3.2390248128063177e-11
        5.7938432761439429e-12
        1.3370330608960321e-12
        3.5262588767213599e-13
        1.0288359890389495e-13
        3.3665560846579745e-14
        1.3003349638950275e-14
        6.1355173976492494e-15
    ]

    for i in eachindex(h)
        ret = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, h[i], F10ₐ; n=n_eq)
        @test ret ≈ results[i] rtol=2e-2
    end
end

@testset "Fetching F10.7 (Modified)" begin
    SpaceIndices.init()

    instant = DateTime("2023-01-01T10:00:00")
    h       = collect(200:200:1000) .* 1000
    ϕ_gd    = deg2rad(-23.0)
    λ       = deg2rad(-45.0)

    # Compute expected with manually supplied F10.7 (81-day avg).
    jd = datetime2julian(instant)
    F10ₐ = sum(space_index.(Val(:F10obs), jd + k) for k in -40:40) / 81

    expected = AtmosphericModels.harrispriester_modified.(instant, ϕ_gd, λ, h, F10ₐ)

    for k in eachindex(h)
        result = AtmosphericModels.harrispriester_modified(instant, ϕ_gd, λ, h[k])
        @test result ≈ expected[k]
    end
end

@testset "Physical Behavior (Modified)" begin
    jd   = date_to_jd(2003, 3, 21, 1, 0, 0)
    ϕ_gd = 0.0
    λ    = 0.0
    n_eq = 2.001

    # Density must decrease with altitude.
    rho_200 = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, 200e3, 150.0; n=n_eq)
    rho_500 = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, 500e3, 150.0; n=n_eq)
    rho_800 = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, 800e3, 150.0; n=n_eq)
    @test rho_200 > rho_500 > rho_800

    # Higher solar activity produces higher density.
    rho_low  = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, 500e3, 75.0;  n=n_eq)
    rho_med  = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, 500e3, 150.0; n=n_eq)
    rho_high = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, 500e3, 250.0; n=n_eq)
    @test rho_low < rho_med < rho_high

    # All densities must be positive.
    for h_km in [120, 200, 400, 600, 800, 1000]
        @test AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, h_km * 1e3, 150.0; n=n_eq) > 0
    end
end

@testset "Near Layer Boundary Continuity (Modified)" begin
    jd   = date_to_jd(2003, 3, 21, 1, 0, 0)
    ϕ_gd = 0.0
    λ    = 0.0
    F10ₐ = 150.0
    n_eq = 2.001

    # FORTRAN reference values near the 300 km layer boundary (±1 km, 0.1 km steps).
    # The Junkins/Jancaitis weighting ensures C³ continuity across boundaries.
    boundary_alts = [
        299.0, 299.1, 299.2, 299.3, 299.4, 299.5, 299.6, 299.7, 299.8, 299.9,
        300.0, 300.1, 300.2, 300.3, 300.4, 300.5, 300.6, 300.7, 300.8, 300.9, 301.0
    ] * 1000

    boundary_rhos = [
        2.2255981257424839e-11
        2.2207493914120607e-11
        2.2159112252435146e-11
        2.2110836041937175e-11
        2.2062665052700210e-11
        2.2014599055297533e-11
        2.1966622283628775e-11
        2.1918639129234863e-11
        2.1870672684999589e-11
        2.1822973961994947e-11
        2.1775835950869826e-11
        2.1729374745502388e-11
        2.1683452074951774e-11
        2.1637776777703074e-11
        2.1592113735725923e-11
        2.1546443980177371e-11
        2.1500856528262700e-11
        2.1455365578194779e-11
        2.1409970925601349e-11
        2.1364672366544803e-11
        2.1319469697517952e-11
    ]

    for i in eachindex(boundary_alts)
        ret = AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, boundary_alts[i], F10ₐ; n=n_eq)
        @test ret ≈ boundary_rhos[i] rtol=2e-2
    end

    # Verify monotonic decrease across the boundary.
    rhos = [AtmosphericModels.harrispriester_modified(jd, ϕ_gd, λ, h, F10ₐ; n=n_eq) for h in boundary_alts]
    @test issorted(rhos; rev=true)
end
