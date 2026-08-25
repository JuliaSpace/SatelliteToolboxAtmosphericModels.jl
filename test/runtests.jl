using Test

using DelimitedFiles
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase

@testset "Atmospheric Model Jacchia 1977" verbose = true begin
    include("./jacchia1977/jacchia1977.jl")
end

@testset "Atmospheric Model Jacchia-Roberts 1971" verbose = true begin
    include("./jr1971.jl")
end

@testset "Atmospheric Model Jacchia-Bowman 2008" verbose = true begin
    include("./jb2008/jb2008.jl")
end

@testset "Atmospheric Model NRLMSISE-00" verbose = true begin
    include("./nrlmsise00.jl")
end

@testset "Exponential Atmospheric Model" verbose = true begin
    include("./exponential.jl")
end

@testset "Harris-Priester Atmospheric Model" verbose = true begin
    include("./harrispriester.jl")
end

if isempty(VERSION.prerelease)
    using Pkg
    Pkg.add("JET")
    Pkg.add("AllocCheck")
    Pkg.add("Aqua")

    using JET
    using AllocCheck
    using Aqua

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end

    # ForwardDiff, Mooncake, ChainRulesCore, and Zygote are already declared in the test
    # target of the Project.toml.
    using ForwardDiff
    using Mooncake
    using ChainRulesCore
    using Zygote

    @testset "Extension Tests" verbose = true begin
        include("./extensions.jl")
    end
else
    @warn "Performance checks and autodiff extensions are not guaranteed to work on julia-nightly, skipping"
end
