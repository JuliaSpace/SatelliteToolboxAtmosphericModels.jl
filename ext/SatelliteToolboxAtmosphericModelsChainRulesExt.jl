module SatelliteToolboxAtmosphericModelsChainRulesExt

using SatelliteToolboxAtmosphericModels

using ChainRulesCore

function ChainRulesCore.rrule(::typeof(AtmosphericModels._get_doy), jd::Number)

    y = AtmosphericModels._get_doy(jd)

    function _get_doy_pullback(Δ::Number)
        return (NoTangent(), Δ)
    end

    return y, _get_doy_pullback

end

function ChainRulesCore.frule((_, Δ), ::typeof(AtmosphericModels._get_doy), jd::Number)

    y = AtmosphericModels._get_doy(jd)

    return y, Δ

end

end