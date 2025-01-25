module SatelliteToolboxAtmosphericModelsMooncakeExt

using SatelliteToolboxAtmosphericModels

using Mooncake
using ChainRulesCore

function ChainRulesCore.rrule(::typeof(AtmosphericModels._get_doy), jd::Number)

    y = AtmosphericModels._get_doy(jd)

    function _get_doy_pullback(Δ::Number)
        return (NoTangent(), Δ)
    end

    return y, _get_doy_pullback

end

Mooncake.@from_rrule DefaultCtx Tuple{typeof(AtmosphericModels._get_doy), Base.IEEEFloat}

end