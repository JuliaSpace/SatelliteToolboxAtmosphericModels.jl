using Test

using DelimitedFiles
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase

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

if isempty(VERSION.prerelease)
    # Add Mooncake and Enzyme to the project if not the nightly version
    # Adding them via the Project.toml isn't working because it tries to compile them before reaching the gating
    using Pkg
    Pkg.add("DifferentiationInterface")
    Pkg.add("Enzyme")
    Pkg.add("FiniteDiff")
    Pkg.add("ForwardDiff")
    Pkg.add("Mooncake")
    Pkg.add("PolyesterForwardDiff")
    Pkg.add("Zygote")

    Pkg.add("JET")
    Pkg.add("AllocCheck")
    Pkg.add("Aqua")

    # Test with Mooncake and Enzyme along with the other backends
    using DifferentiationInterface
    using Enzyme, FiniteDiff, ForwardDiff, Mooncake, PolyesterForwardDiff, Zygote
    const _BACKENDS = (
        ("ForwardDiff", AutoForwardDiff()),
        ("Enzyme", AutoEnzyme()),
        ("Mooncake", AutoMooncake(;config=nothing)),
        ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
        ("Zygote", AutoZygote()),
    )

    using JET
    using AllocCheck
    using Aqua

    @testset "Differentiation Tests" verbose = true begin
        include("./differentiablity.jl")
    end
    
    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end
else
    @warn "Differentiation backends not guaranteed to work on julia-nightly, skipping tests"
    @warn "Performance checks not guaranteed to work on julia-nightly, skipping tests"
end