@testset "JET Testing" begin
    rep = JET.test_package(SatelliteToolboxAtmosphericModels; toplevel_logger=nothing, target_modules=(@__MODULE__,))
end

##########################################################################################################
# These checks don't equal 0 becuase the logger allocates and the flag can't be assessed at compile time
# Checking with BenchmarkTools will show these are actually npn-allocating when verbose=false
# https://julialang.github.io/AllocCheck.jl/dev/tutorials/optional_debugging_and_logging/
##########################################################################################################
rc = zeros(5)
@test length(check_allocs(AtmosphericModels.exponential, (Float64,))) == 0
@test length(check_allocs((x1, x2, x3, x4, rc) -> AtmosphericModels.jr1971(x1, x2, x3, x4; roots_container=rc, verbose=Val(false)), (DateTime, Float64, Float64, Float64, Vector{Float64}))) == 0
@test length(check_allocs((x1, x2, x3, x4) -> AtmosphericModels.jb2008(x1, x2, x3, x4; verbose=Val(false)), (DateTime, Float64, Float64, Float64))) == 0
@test length(check_allocs((x1, x2, x3, x4, P) -> AtmosphericModels.nrlmsise00(x1, x2, x3, x4; P=P, verbose=Val(false)) , (DateTime, Float64, Float64, Float64, Matrix{Float64}))) == 0