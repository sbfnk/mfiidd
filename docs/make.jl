using Documenter
using Documenter.Remotes: GitHub
using MFIIDD

DocMeta.setdocmeta!(MFIIDD, :DocTestSetup, :(using MFIIDD); recursive=true)

makedocs(;
    modules=[MFIIDD],
    sitename="MFIIDD",
    authors="Sebastian Funk and contributors",
    repo=GitHub("sbfnk", "mfiidd"),
    format=Documenter.HTML(;
        canonical="https://sbfnk.github.io/mfiidd/api",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/sbfnk/mfiidd",
    devbranch="main",
    target="build",
    dirname="api",
    push_preview=false,
)
