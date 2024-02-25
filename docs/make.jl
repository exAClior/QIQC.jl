using QIQC
using Documenter

DocMeta.setdocmeta!(QIQC, :DocTestSetup, :(using QIQC); recursive=true)

makedocs(;
    modules=[QIQC],
    authors="Yusheng Zhao <yushengzhao2020@outlook.com> and contributors",
    sitename="QIQC.jl",
    format=Documenter.HTML(;
        canonical="https://exAClior.github.io/QIQC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/exAClior/QIQC.jl",
    devbranch="main",
)
