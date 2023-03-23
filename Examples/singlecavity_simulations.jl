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
