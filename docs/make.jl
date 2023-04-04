push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using WaveguideQED

DocMeta.setdocmeta!(WaveguideQED, :DocTestSetup, :(using WaveguideQED); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"))

makedocs(
bib,
doctest = false,
clean = true,
sitename = "WaveguideQED.jl",
format = Documenter.HTML(assets=["assets/style.css"]),
modules = [WaveguideQED],
authors = "Matias Bundgaard-Nielsen",
pages = [
"WaveguideQED.jl" => "index.md",
"toturial_continous_fockstates.md",
"toturial_combining.md",
"toturial_twochannels.md",
"Beamsplitter interference" => "toturial_detection.md",
"Examples" => ["Input output waveguides" => "example_lodahl.md"],
"API" => "API.md",
"References" => "references.md",
]
)

deploydocs(
    repo = "github.com/mabuni1998/WaveguideQED.jl.git"
)
