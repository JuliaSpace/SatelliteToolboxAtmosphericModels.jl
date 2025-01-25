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

function ChainRulesCore.rrule(::typeof(floor), x::Number)

    y = floor(x)

    function _get_doy_pullback(Δ::Number)
        return (NoTangent(), 0.0)
    end

    return y, _get_doy_pullback

end

Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(floor), Base.IEEEFloat}
Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(AtmosphericModels._get_doy), Base.IEEEFloat}

end