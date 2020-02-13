using Documenter, JellyMe4

makedocs(;
    modules=[JellyMe4],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/palday/JellyMe4.jl/blob/{commit}{path}#L{line}",
    sitename="JellyMe4.jl",
    authors="Phillip Alday",
    assets=String[],
)

deploydocs(;
    repo="github.com/palday/JellyMe4.jl",
)
