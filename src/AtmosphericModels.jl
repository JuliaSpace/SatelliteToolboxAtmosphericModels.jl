## Description #############################################################################
#
# Definition of the module AtmosphericModels to access the models defined here.
#
############################################################################################

module AtmosphericModels

using Accessors
using Crayons
using LinearAlgebra
using PolynomialRoots
using Printf
using SpaceIndices
using SatelliteToolboxBase
using SatelliteToolboxCelestialBodies
using SatelliteToolboxLegendre
using StaticArraysCore

import Base: show

############################################################################################
#                                        Constants                                         #
############################################################################################

const _D = string(Crayon(reset = true))
const _B = string(crayon"bold")

############################################################################################
#                                         Includes                                         #
############################################################################################
include("utils.jl")

include("./exponential/constants.jl")
include("./exponential/exponential.jl")

include("./jr1971/types.jl")
include("./jr1971/constants.jl")
include("./jr1971/jr1971.jl")
include("./jr1971/show.jl")

include("./jb2008/types.jl")
include("./jb2008/constants.jl")
include("./jb2008/jb2008.jl")
include("./jb2008/show.jl")

include("./nrlmsise00/types.jl")
include("./nrlmsise00/auxiliary.jl")
include("./nrlmsise00/constants.jl")
include("./nrlmsise00/math.jl")
include("./nrlmsise00/nrlmsise00.jl")
include("./nrlmsise00/show.jl")

end # module AtmosphericModels
