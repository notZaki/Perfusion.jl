using Documenter, Perfusion

makedocs(
    ;
    modules = [Perfusion],
    format = Documenter.HTML(),
    pages = ["Home" => "index.md",],
    repo = "https://github.com/notZaki/Perfusion.jl/blob/{commit}{path}#L{line}",
    sitename = "Perfusion.jl",
    authors = "Zaki A",
    assets = String[],
)

include("apply_citeproc.jl")
apply_citeproc("./build")

deploydocs(; repo = "github.com/notZaki/Perfusion.jl")
