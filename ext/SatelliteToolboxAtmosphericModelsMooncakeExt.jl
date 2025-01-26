module SatelliteToolboxAtmosphericModelsMooncakeExt

using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase

using Mooncake
using ChainRulesCore

function ChainRulesCore.rrule(::typeof(AtmosphericModels._get_doy), jd::Number)

    y = AtmosphericModels._get_doy(jd)

    function _get_doy_pullback(Δ::Number)
        return (ChainRulesCore.NoTangent(), Δ)
    end

    return y, _get_doy_pullback

end

function ChainRulesCore.rrule(::typeof(jd_to_date), jd::Number)

    y = jd_to_date(jd)

    function jd_to_date_pullback(Δ)
        return (ChainRulesCore.NoTangent(), (0.0, 0.0, 0.0, 0.0, 0.0, Δ))
    end

    return y, jd_to_date_pullback

end

Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(AtmosphericModels._get_doy), Base.IEEEFloat}
Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(jd_to_date), Base.IEEEFloat}

end