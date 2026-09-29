using Test
using WaveguideQED
using QuantumOpticsBase
using QuantumOptics   # waveguide_montecarlo is defined in the QuantumOptics extension

@testset "waveguide_montecarlo" begin
    times = 0:0.1:4
    dt = 0.1
    bw = WaveguideBasis(1, 1, times)
    bc = FockBasis(1)
    H = im * sqrt(1 / dt) * (create(bc) ⊗ destroy(bw) - destroy(bc) ⊗ create(bw))
    ψ0 = fockstate(bc, 1) ⊗ zerophoton(bw)
    n = (create(bc) * destroy(bc)) ⊗ identityoperator(bw)
    fout(t, ψ) = real(expect(n, ψ))
    # without jumps the trajectory follows the Schrödinger equation
    J0 = [0.0 * destroy(bc) ⊗ identityoperator(bw)]
    _, nm = waveguide_montecarlo(times, ψ0, H, J0; fout=fout)
    _, ne = waveguide_evolution(times, ψ0, H; fout=fout)
    @test nm ≈ ne atol = 1e-4
    ψ = waveguide_montecarlo(times, ψ0, H, [destroy(bc) ⊗ identityoperator(bw)]; seed=1)
    @test norm(ψ) ≈ 1
end
