using QuantumOptics
using WaveguideQED
using LinearAlgebra
using PyPlot
pygui(false)
pygui(true)
include("singlecavity_simulations.jl")
#Parameter structure imported from singlecavity_simulations (bit overkill for now, but can prove usefull)
param = parameters()

#Set detuning:
param.δ = 0
#Set non linearity:
param.x3 = 0

param.times = 0:0.1:10

dt = param.times[2] - param.times[1]

#Create operators for two photons interacting with cavity
bc = FockBasis(1)
bw = WaveguideBasis(1,param.times)
a = destroy(bc)
ad = create(bc);
n = ad*a ⊗ identityoperator(bw)
#$w†a and a†w efficient implementation$
wda = emission(bc,bw)
adw = absorption(bc,bw)
H = param.δ*n + im*sqrt(param.γ/dt)*(adw-wda) + param.x3/4*(n*n+n)


#Define input onephoton state shape
ξfun(t,σ,t0) = complex(sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2))
ξvec = ξfun.(param.times,param.σ,param.t0)

ψ_cw = onephoton(bw,ξvec)
psi = fockstate(bc,0) ⊗  ψ_cw 

#Solve
@btime ψ = waveguide_evolution(param.times, psi, H)

#REFERENCE SOLUTION
sol1 = solve_differentialeq(param,ξfun)
ref_sol = ξfun.(sol1.t,param.σ,param.t0)-sqrt(param.γ)*sol1

#Plot single photon waveguide state 
ψ_single = OnePhotonView(ψ)/sqrt(dt)

fig,ax = subplots(1,1,figsize=(9,4.5))
ax.plot(param.times,abs.(ξvec).^2,"g-",label="Input pulse")
ax.plot(param.times,abs.(ψ_single).^2,"ro",label="δ = 0",fillstyle="none")
ax.plot(sol1.t,abs.(ref_sol).^2,"r-")
ax.set_xlabel("time [1/γ]")
ax.set_ylabel(L"$\xi_{out}^{(1)}$")
ax.legend()
plt.tight_layout()

using DifferentialEquations
using LinearAlgebra
using PyPlot
function dpsi!(dpsi,psi,p,t)
    y,d,nsteps,dt = p
    timeindex = round(Int,t/dt,RoundDown) + 1
    dpsi .= 0
    dpsi[2+nsteps] = -sqrt(y/dt)*psi[1+timeindex]
    dpsi[1+timeindex] = sqrt(y/dt)*psi[2+nsteps]
end

#Define parameters
y,d,dt = 1,0,0.1
times = 0:dt:10
N = length(times)
p = (y,d,N,dt)

#Define input gaussian state with width s = 1 and arrivial time t0=5
xi(t,s,t0) = complex(sqrt(2/s)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/s^2))
psi = zeros(ComplexF64,2*(N+1))
psi[2:N+1] .= sqrt(dt)*xi.(times,1,5)
prob = ODEProblem(dpsi!, psi, (times[1], times[end]+dt), p)
sol = solve(prob, OrdinaryDiffEq.DP5();reltol = 1.0e-8,abstol = 1.0e-10);
psi_out = sol[end][2:N+1]

#Plot the output
fig,ax = subplots(1,1,figsize=(9,4.5))
ax.plot(times,abs.(psi_out).^2,"r-",label="Input pulse")