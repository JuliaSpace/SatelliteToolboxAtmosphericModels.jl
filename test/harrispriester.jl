## Description #############################################################################
#
# Tests related to the Harris-Priester models.
#
############################################################################################

# == Functions: harrispriester =============================================================

@testset "Harris-Priester Model" begin
    # Common inputs
    jd      = date_to_jd(2003, 3, 21, 1, 0, 0)
    instant = julian2datetime(jd)
    
    # Position at 500 km height above the Equator
    ϕ_gd    = 0.0
    λ       = 0.0
    h       = 500e3 # 500 km

    # == Standard Test =====================================================================
    @testset "Standard" begin
        # Expected density from Java test: 3.9237E-13 kg/m^3
        expected_rho = 3.9237e-13
        
        rho = AtmosphericModels.harrispriester(instant, ϕ_gd, λ, h)
        @test rho ≈ expected_rho atol=1e-17
        
        rho_jd = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h)
        @test rho_jd ≈ expected_rho atol=1e-17
    end

    # == Parameter N Test ==================================================================
    @testset "Parameter N" begin
        # Java test logic:
        # rho4 = density with n=4
        # rho2 = density with n=2
        # rho6 = density with n=6
        # c2Psi2 = 0.02163787
        # Assert: (rho6-rho2)/(rho4-rho2) - 1 approx c2Psi2
        
        rho4 = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h; n=4)
        rho2 = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h; n=2)
        rho6 = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h; n=6)
        
        c2Psi2 = 0.02163787
        val = (rho6 - rho2) / (rho4 - rho2) - 1.0
        # Note: Small discrepancy (~0.68%) due to differences in Sun position/GMST calculations
        # between Orekit (Java) and SatelliteToolbox (Julia)
        @test val ≈ c2Psi2 rtol=1e-2
    end

    # == Max Altitude Test =================================================================
    @testset "Max Altitude" begin
        # Position at 1500 km height
        h_high = 1500e3
        
        rho = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h_high)
        @test rho == 0.0
    end

    # == User Table Test ===================================================================
    @testset "User Table" begin
        # Custom table from Java test
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

        h_high = 1500e3
        
        # Note: the Java test uses `earth.transform(point)` which handles coordinate
        # transformation. Our `harrispriester` takes `ϕ_gd, λ, h`.
        # For the user table test, it seems to calculate density at 1500 km.
        # Expected: 2.9049E-7 kg/m^3
        expected_rho = 2.9049e-7
        
        rho = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h_high; alt_ρ=user_tab)
        @test rho ≈ expected_rho atol=1e-11

        rho6 = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h_high; n=6, alt_ρ=user_tab)
        rho2 = AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h_high; n=2, alt_ρ=user_tab)
        
        c2Psi2 = 0.02163787
        val = (rho6 - rho2) / (rho - rho2) - 1.0
        # Note: Small discrepancy (~0.68%) due to differences in Sun position/GMST calculations
        # between Orekit (Java) and SatelliteToolbox (Julia)
        @test val ≈ c2Psi2 rtol=1e-2
    end

    # == Out of Range Test =================================================================
    @testset "Out of Range" begin
        h_low = 50e3 # 50 km
        @test_throws ArgumentError AtmosphericModels.harrispriester(jd, ϕ_gd, λ, h_low)
    end
end


@testset "Harris-Priester Model - Modified (Hatten & Russell 2017)" begin
    # Tests comparing Julia modified Harris-Priester against FORTRAN reference
    # from harris_priester_mod_dist.f90 (Hatten & Russell 2017).
    #
    # Test conditions: JD=2452720.541666667 (2003-03-21 01:00:00 UTC), F10.7=150 SFU
    #
    # IMPORTANT NOTES:
    # 1. The FORTRAN code computes n (npow) from orbital inclination:
    #    npow = 2.001 + 4*sin²(inclination)
    # 2. FORTRAN uses ECI position, Julia uses geodetic coordinates.
    #    Geodetic longitude must be calculated from ECI position at the test epoch.
    # 3. For equatorial orbits (i=0°): n = 2.001
    
    # Common test parameters
    jd = date_to_jd(2003, 3, 21, 1, 0, 0)
    instant = julian2datetime(jd)
    F10_avg = 150.0  # 81-day average
    
    # n value for equatorial orbit (from FORTRAN: npow = 2.001 + 4*sin²(0°) = 2.001)
    n_equatorial = 2.001
    
    @testset "500 km, Equator" begin
        h = 500e3
        ϕ_gd = 0.0
        λ = deg2rad(165.87)
        
        # FORTRAN reference value
        fortran_rho = 1.367889737809932e-12  # kg/m³
        
        # Julia implementation with correct n for equatorial orbit
        julia_rho = AtmosphericModels.harrispriester_modified(instant, ϕ_gd, λ, h, F10_avg; n=n_equatorial)
        
        # Test against FORTRAN reference
        @test julia_rho ≈ fortran_rho rtol=1e-2
    end
    
    @testset "800 km, Equator" begin
        h = 800e3
        ϕ_gd = 0.0
        λ = deg2rad(165.87)
        
        # FORTRAN reference value
        fortran_rho = 3.454417802095930e-14  # kg/m³
        
        julia_rho = AtmosphericModels.harrispriester_modified(instant, ϕ_gd, λ, h, F10_avg; n=n_equatorial)
        
        @test julia_rho ≈ fortran_rho rtol=1e-2
    end
    
    @testset "200 km, Equator" begin
        h = 200e3
        ϕ_gd = 0.0
        λ = deg2rad(165.87)
        
        # FORTRAN reference value
        fortran_rho = 3.216921952487804e-10  # kg/m³
        
        julia_rho = AtmosphericModels.harrispriester_modified(instant, ϕ_gd, λ, h, F10_avg; n=n_equatorial)
        
        @test julia_rho ≈ fortran_rho rtol=1e-2
    end
    
    @testset "Physical behavior" begin
        λ = deg2rad(165.87)
        
        # Density must decrease with altitude
        rho_200 = AtmosphericModels.harrispriester_modified(instant, 0.0, λ, 200e3, F10_avg; n=n_equatorial)
        rho_500 = AtmosphericModels.harrispriester_modified(instant, 0.0, λ, 500e3, F10_avg; n=n_equatorial)
        rho_800 = AtmosphericModels.harrispriester_modified(instant, 0.0, λ, 800e3, F10_avg; n=n_equatorial)
        @test rho_200 > rho_500 > rho_800
    end
end