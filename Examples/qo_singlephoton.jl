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
param.δ = 2
#Set non linearity:
param.x3 = 0

param.times = 0:0.1:20

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
ψ = waveguide_evolution(param.times, psi, H)

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

