using CoinfectionSimulator
using Documenter

DocMeta.setdocmeta!(CoinfectionSimulator, :DocTestSetup, :(using CoinfectionSimulator); recursive=true)

makedocs(;
    modules=[CoinfectionSimulator],
    authors="July Pilowsky <japilo@users.noreply.github.com> and contributors",
    repo="https://github.com/japilo/CoinfectionSimulator.jl/blob/{commit}{path}#{line}",
    sitename="CoinfectionSimulator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://japilo.github.io/CoinfectionSimulator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
    ],
)

deploydocs(;
    repo="github.com/japilo/CoinfectionSimulator.jl",
    devbranch="main",
)
