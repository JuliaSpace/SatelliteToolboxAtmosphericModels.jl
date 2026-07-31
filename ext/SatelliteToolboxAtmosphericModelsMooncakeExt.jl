## Description #############################################################################
#
# This file integrates the atmospheric models with Mooncake.jl. Notice that the rule for
# `AtmosphericModels._get_doy` is defined in the ChainRulesCore.jl extension, which is
# always loaded together with this one.
#
############################################################################################

module SatelliteToolboxAtmosphericModelsMooncakeExt

using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase

using Mooncake
using ChainRulesCore

Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(AtmosphericModels._get_doy), Number}
Mooncake.@zero_adjoint Mooncake.DefaultCtx Tuple{typeof(jd_to_date), Number}

end
