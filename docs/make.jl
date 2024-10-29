using ProjMPSs
using Documenter

DocMeta.setdocmeta!(ProjMPSs, :DocTestSetup, :(using ProjMPSs); recursive=true)

makedocs(;
    modules=[ProjMPSs],
    authors="Hiroshi Shinaoka <h.shinaoka@gmail.com> and contributors",
    sitename="ProjMPSs.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/tensor4all/ProjMPSs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/tensor4all/ProjMPSs.jl.git", devbranch="main")
