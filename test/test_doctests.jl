using Documenter
using WaveguideQED

@testset "Doctests" begin
    DocMeta.setdocmeta!(WaveguideQED, :DocTestSetup, :(using WaveguideQED); recursive=true)
    doctest(WaveguideQED)
end
