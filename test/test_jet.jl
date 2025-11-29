@testset "JET checks" begin
    using JET
    using Test
    using WaveguideQED

    rep = JET.report_package(WaveguideQED, target_defined_modules = true)
    println(rep)
    @testset length(JET.get_reports(rep)) <= 10
    @test_broken length(JET.get_reports(rep)) == 0
end