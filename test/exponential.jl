## Description #############################################################################
#
# Tests related to the exponential atmospheric model.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
############################################################################################

# == Functions: exponential ================================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# == Scenario 01 ===========================================================================
#
#   Example 8-4 [1, p. 567]
#
#   Using the exponential atmospheric model, one gets:
#
#       ρ(747.2119 km) = 2.1219854 ⋅ 10⁻⁴ kg / m³
#
# == Scenario 02 ===========================================================================
#
#   By definition, one gets `ρ(h₀) = ρ₀(h₀)` for all tabulated values.
#
# == Scenario 03 ===========================================================================
#
#   Inside every interval, the atmospheric density must be monotonically
#   decreasing.
#
############################################################################################

@testset "Default Tests" begin
    # == Scenario 01 =======================================================================

    @test AtmosphericModels.exponential(747211.9) ≈ 2.1219854e-14 rtol = 1e-8

    # == Scenario 02 =======================================================================

    for i in 1:length(AtmosphericModels._EXPONENTIAL_ATMOSPHERE_H₀)
        h = 1000 * AtmosphericModels._EXPONENTIAL_ATMOSPHERE_H₀[i]

        @test AtmosphericModels.exponential(h) == AtmosphericModels._EXPONENTIAL_ATMOSPHERE_ρ₀[i]
    end

    # == Scenario 03 =======================================================================

    for i in 2:length(AtmosphericModels._EXPONENTIAL_ATMOSPHERE_H₀)
        h₀i = 1000 * AtmosphericModels._EXPONENTIAL_ATMOSPHERE_H₀[i - 1]
        h₀f = 1000 * AtmosphericModels._EXPONENTIAL_ATMOSPHERE_H₀[i]
        Δ   = 100
        Δh₀ = (h₀f - h₀i)/Δ

        ρk = AtmosphericModels.exponential(h₀i)

        for k in 2:Δ
            ρk₋₁ = ρk
            h    = h₀i + Δh₀ * (k - 1)
            ρk   = AtmosphericModels.exponential(h)

            @test ρk < ρk₋₁
        end
    end
end
