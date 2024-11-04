using JET
#@testset "JET Testing" begin
    rep = JET.test_package(SatelliteToolboxAtmosphericModels; toplevel_logger=nothing)
#end

reps = report_package(SatelliteToolboxAtmosphericModels)

report = JET.get_reports(reps)

for i in 1:19
    @report_opt target_modules=(@__MODULE__,) AtmosphericModels.jr1971(instant, ϕ_gd, λ, h[i])
end

using Cthulhu

report[1]

using AllocCheck
using SatelliteToolboxAtmosphericModels

instant = DateTime("2023-01-01T10:00:00")
h       = collect(90:50:1000) .* 1000
ϕ_gd    = -23 |> deg2rad
λ       = -45 |> deg2rad
F10     = 152.6
F10ₐ    = 159.12345679012347
Kp      = 2.667

using SpaceIndices
SpaceIndices.init()

using Test
using BenchmarkTools

##########################################################################################################
# These checks don't equal 0 becuase the logger allocates and the flag can't be assessed at compile time
# Checking with BenchmarkTools will show these are actually npn-allocating when verbose=false
# https://julialang.github.io/AllocCheck.jl/dev/tutorials/optional_debugging_and_logging/
##########################################################################################################
@test length(check_allocs(AtmosphericModels.exponential, (Float64,))) == 0
@test length(check_allocs(AtmosphericModels.jr1971, (DateTime, Float64, Float64, Float64))) == 133
@test length(check_allocs(AtmosphericModels.jb2008, (DateTime, Float64, Float64, Float64))) == 135
@test length(check_allocs((x1, x2, x3, x4, P) -> AtmosphericModels.nrlmsise00(x1, x2, x3, x4; P=P) , (DateTime, Float64, Float64, Float64, Matrix{Float64}))) == 138
