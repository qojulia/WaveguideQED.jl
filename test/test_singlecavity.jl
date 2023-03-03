include("helper_functions.jl")

@testset "Singlephoton with cavity" begin
    #Set parameters
    param = parameters()
    param.δ = 0
    param.x3 = 0
    param.times = 0:0.1:20
    dt = param.times[2] - param.times[1]

    #Create operators for two photons interacting with cavity
    bc = FockBasis(1)
    bw = WaveguideBasis(1,param.times)
    btotal = tensor(bc,bw)
    a = destroy(bc)
    ad = create(bc);
    n = ad*a ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w
    H = param.δ*n + im*sqrt(param.γ/dt)*(adw-wda) + param.x3/4*(n*n+n)
    #Define input twophoton state shape
    ξfun(t,σ,t0) = complex( sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2))
    ξvec = ξfun.(param.times,param.σ,param.t0)
    #Define initial state
    ψ_cw = onephoton(bw,ξvec)
    psi = fockstate(bc,0) ⊗  ψ_cw 

    #Solve
    ψ = waveguide_evolution(param.times, psi, H)

    #REFERENCE SOLUTION
    sol1 = solve_differentialeq(param,ξfun)
    ref_sol = ξfun.(sol1.t,param.σ,param.t0)-sqrt(param.γ)*sol1

    ψ_single = view_onephoton(ψ)/sqrt(dt)
    @test isapprox(ψ_single,ref_sol,rtol=0.05)
end