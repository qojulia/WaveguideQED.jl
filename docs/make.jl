push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using WaveguideQED

DocMeta.setdocmeta!(WaveguideQED, :DocTestSetup, :(using WaveguideQED); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"))

makedocs(
bib,
doctest = true,
clean = true,
sitename = "WaveguideQED.jl",
format = Documenter.HTML(assets=["assets/style.css"]),
modules = [WaveguideQED],
authors = "Matias Bundgaard-Nielsen",
pages = [
"WaveguideQED.jl" => "index.md",
"theoreticalbackground.md",
"tutorial.md",
"multiplewaveguides_v2.md",
"beamsplitter.md",
"input_output.md",
"time_delay.md",
#"Examples" => ["Scattering on two level system" => "example_lodahl.md"],
"API" => "API.md",
"References and suggested readings" => "references.md",
"Experimental or Outdated functionalities" => ["detection.md"]
]
)

deploydocs(
    repo = "github.com/qojulia/WaveguideQED.jl.git"
)
