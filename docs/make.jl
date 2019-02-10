using Documenter, TightBinding

makedocs(;
    modules=[TightBinding],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/cometscome/TightBinding.jl/blob/{commit}{path}#L{line}",
    sitename="TightBinding.jl",
    authors="Yuki Nagai",
    assets=[],
)

deploydocs(;
    repo="github.com/cometscome/TightBinding.jl",
)
