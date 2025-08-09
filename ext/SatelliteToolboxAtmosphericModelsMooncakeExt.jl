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

Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(AtmosphericModels._get_doy), Number}
Mooncake.@zero_adjoint Mooncake.DefaultCtx Tuple{typeof(jd_to_date), Number}

end