using ImpactChron
using Documenter

DocMeta.setdocmeta!(ImpactChron, :DocTestSetup, :(using ImpactChron); recursive=true)

makedocs(;
    modules=[ImpactChron],
    authors="Graham Harper Edwards",
    repo="https://github.com/grahamedwards/ImpactChron.jl/blob/{commit}{path}#{line}",
    sitename="ImpactChron.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://grahamedwards.github.io/ImpactChron.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/grahamedwards/ImpactChron.jl",
    devbranch="main",
)
