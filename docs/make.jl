using sfED
using Documenter

push!(:OAD_PATH, "../src/")
makedocs(;
    modules=[sfED],
    authors=["Steffen Backes <mail@mail.mm> and Contributers"],
    repo="https://github.com/steffenbackes/sfED/blob/{commit}{path}@L{line}",
    sitename="TODO: details here",
    format=Documenter.HTML(;
        prettyurl=get(ENV, "CI", nothing) == "true",
        canonical="https://github.com/steffenbackes/sfED",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    branch="hg-pages",
    devbranch="main",
    devurl ="stable",
    repo="https://github.com/steffenbackes/sfED"
)
