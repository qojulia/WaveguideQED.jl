using LinearAlgebra
include("helper_functions.jl")

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
    testvec = ones(length(param.times)) .* 1/sqrt(get_nsteps(psi.basis)*(get_nsteps(psi.basis)-1)/2+sqrt(2)*get_nsteps(psi.basis))
    @test isapprox(ψ_single,testvec,rtol=0.01)
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
    testvec[get_waveguidetimeindex(wda)] *= 2
    ψ_double=TwophotonTimestepView(view_waveguide(psi_out),get_waveguidetimeindex(wda),get_nsteps(psi.basis))
    @test isapprox([ψ_double[i] for i in eachindex(ψ_double)],testvec)    
end

@testset "1 photon 2 photon the same" begin
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
    psi = fockstate(bc,0) ⊗  input
    wda_out_2 = copy(psi)
    adw_out_2 = copy(psi)
    QuantumOpticsBase.mul!(wda_out_2,wda,psi,1,0.0)
    QuantumOpticsBase.mul!(adw_out_2,wda,psi,1,0.0)

    bw = WaveguideBasis(1,param.times)
    bc = FockBasis(2)
    a = destroy(bc)
    ad = create(bc);
    n = ad*a ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w
    
    input = onephoton(bw,ξfun,param.times)
    psi = fockstate(bc,0) ⊗  input
    adw_out_1 = copy(psi)
    wda_out_1 = copy(psi)
    QuantumOpticsBase.mul!(wda_out_1,wda,psi,1,0.0)
    QuantumOpticsBase.mul!(adw_out_1,adw,psi,1,0.0)
    
    @test isapprox(view_onephoton(wda_out_1),view_onephoton(wda_out_2))
    @test isapprox(view_onephoton(adw_out_1),view_onephoton(adw_out_2))  
    @test isapprox(wda_out_1.data[1],wda_out_2.data[1])
    @test isapprox(adw_out_1.data[1],adw_out_2.data[1]) 
end

@testset "Cavity Waveguide Operator" begin
    param = parameters()
    param.δ = 0
    param.x3 = 0

    param.times = 0:0.1:10
    dt = param.times[2] -param.times[1]
    bw = WaveguideBasis(2,param.times)
    bc = FockBasis(1)
    be = FockBasis(1)
    a = destroy(bc)
    ad = create(bc);
    n_l = identityoperator(be) ⊗ (ad*a) ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);

    wda_l = identityoperator(be) ⊗ emission(bc,bw)
    adw_l = identityoperator(be) ⊗ absorption(bc,bw)
    H_l = param.δ*n_l+im*sqrt(param.γ/dt)*(adw_l-wda_l)+param.x3/4*n_l*n_l-param.x3/4*n_l


    n_f = identityoperator(bw) ⊗ (ad*a) ⊗ identityoperator(be)
    wda_f = emission(bw,bc) ⊗ identityoperator(be)
    adw_f = absorption(bw,bc) ⊗ identityoperator(be)
    H_f = param.δ*n_f+im*sqrt(param.γ/dt)*(adw_f-wda_f)+param.x3/4*n_f*n_f-param.x3/4*n_f


    ξfun(t,σ,t0) = complex(1/(σ*sqrt(2*pi))*exp(-2*log(2)*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
    ξvec = ξfun.(param.times,1,1) * transpose(ξfun.(param.times,1,1))
    #Define initial state
    #ψ_cw = onephoton(bw,ξfun,param.times,param.σ,param.t0)
    ψ_cw = twophoton(bw,ξvec)

    psi_l = fockstate(be,0) ⊗ fockstate(bc,0) ⊗ ψ_cw
    out_l_emission = copy(psi_l)
    out_l_wda = copy(psi_l)
    out_l_absorption = copy(psi_l)
    out_l_adw = copy(psi_l)


    psi_f = ψ_cw ⊗ fockstate(bc,0) ⊗ fockstate(be,0)
    out_f_emission = copy(psi_f)
    out_f_wda = copy(psi_f)
    out_f_absorption = copy(psi_f)
    out_f_adw = copy(psi_f)


    set_waveguidetimeindex!([wda_l,wda_f,adw_l,adw_f],100)

    mul!(out_l_emission,wda_l,psi_l,1.0,0.0)
    mul!(out_f_emission,wda_f,psi_f,1.0,0.0)
    mul!(out_l_absorption,adw_l,psi_l,1.0,0.0)
    mul!(out_f_absorption,adw_f,psi_f,1.0,0.0)

    wda_l = identityoperator(be) ⊗ (a ⊗ wd)
    adw_l = identityoperator(be) ⊗ (ad ⊗ w)
    H_l = param.δ*n_l+im*sqrt(param.γ/dt)*(adw_l-wda_l)+param.x3/4*n_l*n_l-param.x3/4*n_l

    n_f = identityoperator(bw) ⊗ (ad*a) ⊗ identityoperator(be)
    wda_f = wd ⊗ a ⊗ identityoperator(be)
    adw_f = w ⊗ ad ⊗ identityoperator(be)
    H_f = param.δ*n_f+im*sqrt(param.γ/dt)*(adw_f-wda_f)+param.x3/4*n_f*n_f-param.x3/4*n_f

    set_waveguidetimeindex!([wda_l,wda_f,adw_l,adw_f],100)

    mul!(out_l_wda,wda_l,psi_l,1.0,0.0)
    mul!(out_f_wda,wda_f,psi_f,1.0,0.0)
    mul!(out_l_adw,adw_l,psi_l,1.0,0.0)
    mul!(out_f_adw,adw_f,psi_f,1.0,0.0)

    @test isapprox(out_l_wda,out_l_emission)
    @test isapprox(out_l_adw,out_l_absorption)
    @test isapprox(out_f_wda,out_f_emission)
    @test isapprox(out_f_adw,out_f_absorption)
end