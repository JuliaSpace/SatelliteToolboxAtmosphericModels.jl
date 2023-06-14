using Documenter
using SatelliteToolboxAtmosphericModels

makedocs(
    modules = [SatelliteToolboxAtmosphericModels],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteToolboxAtmosphericModels.jl/stable/",
    ),
    sitename = "Satellite Toolbox Atmospheric Models",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Atmospheric Models" => [
            "Exponential" => "man/exponential.md",
            "Jacchia-Roberts 1971" => "man/jr1971.md",
            "Jacchia-Bowman 2008" => "man/jb2008.md",
            "NRLMSISE-00" => "man/nrlmsise00.md"
        ],
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl.git",
    target = "build",
)
