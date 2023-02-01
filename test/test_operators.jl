using LinearAlgebra
include("paramfile.jl")

@testset "Waveguide annihilation" begin
    param=parameters()
    param = parameters()
    param.δ = 0
    param.x3 = 0
    param.times = 0:0.1:20
    dt = param.times[2] - param.times[1]

    bc = FockBasis(2)
    bw = WaveguideBasis(2,param.times)
    a = destroy(bc)
    ad = create(bc);
    n = ad*a ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w
    
    ξfun(t1,t2) = 1
    input = twophoton(bw,ξfun,param.times)
    psi = fockstate(bc,0) ⊗ input
    psi_out = copy(psi)

    QuantumOpticsBase.mul!(psi_out,adw,psi,1.0,0.0)
    ψ_single=view_onephoton(psi_out,[2,:])
    testvec = ones(length(param.times)) .* 1/sqrt(get_nsteps(psi.basis)*(get_nsteps(psi.basis)+1)/2)
    testvec[get_waveguidetimeindex(psi.basis)] *= sqrt(2)
    @test isapprox(ψ_single,testvec)
end

@testset "Waveguide creation" begin
    param=parameters()
    param = parameters()
    param.δ = 0
    param.x3 = 0
    param.times = 0:0.1:20
    dt = param.times[2] - param.times[1]

    bw = WaveguideBasis(2,param.times)
    bc = FockBasis(2)
    a = destroy(bc)
    ad = create(bc);
    n = ad*a ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w
    
    ξfun(t1) = 1
    input = onephoton(bw,ξfun,param.times)
    psi = fockstate(bc,1) ⊗  input
    psi_out = copy(psi)
    QuantumOpticsBase.mul!(psi_out,wda,psi,1,0.0)
    testvec = ones(length(param.times)) .* 1/sqrt(get_nsteps(psi.basis))
    testvec[get_waveguidetimeindex(psi.basis)] *= 2
    ψ_double=TwophotonView(view_waveguide(psi_out),get_waveguidetimeindex(psi_out.basis),get_nsteps(psi.basis))
    @test isapprox([ψ_double[i] for i in eachindex(ψ_double)],testvec)    
end
