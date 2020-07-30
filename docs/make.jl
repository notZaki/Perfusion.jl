using Documenter, Perfusion

makedocs(
    ;
    modules = [Perfusion],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages=[
        "Introduction" => "index.md",
        "Pre-processing" => Any[
            "Loading DICOM" => "preprocessing/loading.md",
            "T1 mapping" => "preprocessing/T1mapping.md",
            "Signal to Concentration" => "preprocessing/signal2concentration.md"
        ],
        "AIF-based Models" => Any[
            "AIF" => "models/aif.md",
            "Tofts" => "models/tofts.md",
            "Extended Tofts" => "models/extendedtofts.md",
            "Exchange" => "models/exchange.md",
            "Uptake" => "models/uptake.md",
            "Filtration" => "models/filtration.md"
        ],
        "Functions" => "functions.md"
    ],
    repo = "https://github.com/notZaki/Perfusion.jl/blob/{commit}{path}#L{line}",
    sitename = "Perfusion.jl",
)

# Apply citations
include("apply_citeproc.jl")
apply_citeproc("./build")

deploydocs(; repo = "github.com/notZaki/Perfusion.jl")
