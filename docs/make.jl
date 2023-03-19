using FusionSpecies
using Documenter

DocMeta.setdocmeta!(FusionSpecies, :DocTestSetup, :(using FusionSpecies); recursive=true)

makedocs(;
    modules=[FusionSpecies],
    authors="Jerome Guterl",
    repo="https://github.com/jguterl/FusionSpecies.jl/blob/{commit}{path}#{line}",
    sitename="FusionSpecies.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jguterl.github.io/FusionSpecies.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jguterl/FusionSpecies.jl",
    devbranch="main",
)
