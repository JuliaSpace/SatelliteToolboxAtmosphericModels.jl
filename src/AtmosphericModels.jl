## Description #############################################################################
#
# Definition of the module AtmosphericModels to access the models defined here.
#
############################################################################################

module AtmosphericModels

using Dates
using Printf

using Accessors
using SatelliteToolboxBase
using SatelliteToolboxCelestialBodies
using SatelliteToolboxLegendre
using SpaceIndices
using StaticArraysCore

import Base: show

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./types.jl")
include("./utils.jl")

include("./exponential/constants.jl")
include("./exponential/exponential.jl")

include("./harrispriester/constants.jl")
include("./harrispriester/harrispriester.jl")
include("./harrispriester/harrispriester_modified.jl")

include("./jacchia1977/types.jl")
include("./jacchia1977/constants.jl")
include("./jacchia1977/jacchia1977.jl")

include("./jr1971/types.jl")
include("./jr1971/constants.jl")
include("./jr1971/jr1971.jl")

include("./jb2008/types.jl")
include("./jb2008/constants.jl")
include("./jb2008/jb2008.jl")

include("./nrlmsise00/types.jl")
include("./nrlmsise00/auxiliary.jl")
include("./nrlmsise00/constants.jl")
include("./nrlmsise00/math.jl")
include("./nrlmsise00/nrlmsise00.jl")

include("./show.jl")

end # module AtmosphericModels
