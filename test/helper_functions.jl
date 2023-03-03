using Parameters
using DifferentialEquations

@with_kw mutable struct parameters
    γ::Float64 = 1.0
    δ::Float64 = 0.0
    σ::Float64 = 1.0
    t0::Float64 = 5.0
    times = 0.0:0.1:20.0
    x3::Float64 = 0.0
    δb::Float64 = 0.0
    Λn::Float64 = 0.0
end

function solve_differentialeq(param::parameters,ξfun)
    @unpack γ,δ,σ,t0,times = param
    p = [δ,γ]
    f(ψ,p,t) = complex(-(im*p[1]+p[2]/2)*ψ + sqrt(p[2])*ξfun(t,σ,t0))
    prob = ODEProblem(f,sqrt(p[2])*ξfun(times[1],σ,t0),(times[1],times[end]),p)
    solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,saveat=times)
end


@with_kw mutable struct BarretKokParameters
    γ::Float64 = 0.01
    g1::Float64 = 0.3
    g2::Float64 = 0.3
    κ1::Float64 = 1
    κ2::Float64 = 1
    δ::Float64 = 0.0
    td::Float64 = 5.0
    twait::Float64 = 30.0
    trelax::Float64 = 30.0
    times = 0.0:0.5:30.0
    montecarlo_repitions::Int = 10
end

function prep_fast(param)
    param.times = 0:0.5:param.twait
    dt = param.times[2] - param.times[1]
    param.td = param.twait

    #Create operators for two photons interacting with cavity
    bc = FockBasis(1)
    be = NLevelBasis(3)
    bw = WaveguideBasis(1,param.times)
    b_a = bw ⊗ bc ⊗ be
    a = embed(b_a,2,destroy(bc))
    ad = embed(b_a,2,create(bc))
    e = embed(b_a,3,transition(be,1,3))
    ed = embed(b_a,3,transition(be,3,1))
    wda = emission(bw,bc) ⊗ identityoperator(be)
    adw = absorption(bw,bc) ⊗ identityoperator(be)

    H_a = im*sqrt(param.κ1/2/dt)*(adw-wda) + param.g1/2*(ad*e + a*ed) 
    H_b = im*sqrt(param.κ2/2/dt)*(adw-wda)  + param.g2/2*(ad*e + a*ed)

    jump_a = sqrt(param.γ) *identityoperator(bw) ⊗ identityoperator(bc) ⊗ transition(be,1,3)
    jump_b = sqrt(param.γ) *identityoperator(bw) ⊗ identityoperator(bc) ⊗ transition(be,1,3) 

    wa = LazyTensor(b_a,1,destroy(bw))
    wb = LazyTensor(b_a,1,destroy(bw))

    n_w = LazyProduct(dagger(wa),wa)
    
    proj_up = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,2)
    proj_down = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,1)
    p1 = LazyTensorKet(proj_up,proj_down)
    p2 = LazyTensorKet(proj_down,proj_up)
    projector_plus = (p1+p2)/sqrt(2)
    projector_minus = (p1-p2)/sqrt(2)
  
    psi_a = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ (1/sqrt(2)*(nlevelstate(be,2)+nlevelstate(be,3)))
    psi_b = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ (1/sqrt(2)*(nlevelstate(be,2)+nlevelstate(be,3)))
    


    return (;(Base.@locals)...)
end
