using sfED
using Documenter

push!(LOAD_PATH, "../src/")
makedocs(;
    modules=[sfED],
    authors="Steffen Backes <mail@mail.mm> and Contributers",
    repo="https://github.com/steffenbackes/sfED/blob/{commit}{path}#L{line}",
    sitename="TODO: details here",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://steffenbackes.github.io/sfED",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    branch="gh-pages",
    devbranch="main",
    devurl ="stable",
    repo="github.com/steffenbackes/sfED.git"
)
