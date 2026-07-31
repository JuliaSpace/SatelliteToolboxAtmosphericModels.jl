## Description #############################################################################
#
# This file defines ChainRulesCore.jl rules shared by all reverse-mode automatic
# differentiation backends (e.g. Zygote.jl and Mooncake.jl).
#
############################################################################################

module SatelliteToolboxAtmosphericModelsChainRulesCoreExt

using SatelliteToolboxAtmosphericModels

using ChainRulesCore

# The function `_get_doy` computes the day of the year using operations that are not
# differentiable (e.g. `Dates` arithmetic). However, the day of the year is an affine
# function of the Julian day with unitary slope inside a year. Hence, we define its
# derivative as 1.
function ChainRulesCore.rrule(::typeof(AtmosphericModels._get_doy), jd::Number)
    y = AtmosphericModels._get_doy(jd)

    function _get_doy_pullback(Δ)
        return (NoTangent(), Δ)
    end

    return y, _get_doy_pullback
end

end
