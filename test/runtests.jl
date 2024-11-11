using Test

using DelimitedFiles
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase

using DifferentiationInterface
using FiniteDiff, ForwardDiff, Diffractor, Enzyme, Mooncake, PolyesterForwardDiff, ReverseDiff, Tracker, Zygote

using AllocCheck
using JET

@testset "Atmospheric Model Jacchia-Roberts 1971" verbose = true begin
    include("./jr1971.jl")
end

@testset "Atmospheric Model Jacchia-Bowman 2008" verbose = true begin
    cd("./jb2008")
    include("./jb2008/jb2008.jl")
    cd("..")
end

@testset "Atmospheric Model NRLMSISE-00" verbose = true begin
    include("./nrlmsise00.jl")
end

@testset "Exponential Atmospheric Model" verbose = true begin
    include("./exponential.jl")
end

@testset "Differentiation Tests" verbose = true begin
    include("./differentiablity.jl")
end

@testset "Performance Tests" verbose = true begin
    include("./allocations.jl")
end
