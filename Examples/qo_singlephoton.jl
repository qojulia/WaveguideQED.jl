using QuantumOptics
using CavityWaveguide
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

param.times = 0:0.1:20

dt = param.times[2] - param.times[1]
tend = param.times[end]

#Create operators for two photons interacting with cavity
bc = FockBasis(1)
bw = WaveguideBasis(1,param.times)
btotal = tensor(bc,bw)
a = sparse(destroy(bc))
ad = sparse(create(bc));
n = ad*a ⊗ identityoperator(bw)
w = destroy(bw)
wd = create(bw);
wda = a ⊗ wd
adw = ad ⊗ w
H = param.δ*n + im*sqrt(param.γ/dt)*(adw-wda) + param.x3/4*(n*n+n)


#Define input twophoton state shape
ξfun(t,σ,t0) = complex(1/(σ*sqrt(2*pi))*exp(-2*log(2)*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
ξvec = ξfun.(param.times,param.σ,param.t0)
#Define initial state
#ψ_cw = onephoton(bw,ξfun,param.times,param.σ,param.t0)
ψ_cw = onephoton(bw,ξvec)
psi = fockstate(bc,0) ⊗  ψ_cw 

#Solve
ψ = waveguide_evolution(param.times, psi, H)

#REFERENCE SOLUTION
solve_differentialeq(param,ξfun)
sol1 = solve_differentialeq(param,ξfun)
ref_sol = ξfun.(sol1.t,param.σ,param.t0)-sqrt(param.γ)*sol1

#Plot single photon waveguide state 
ψ_single = view_onephoton(ψ)/sqrt(dt)

fig,ax = subplots(1,1,figsize=(9,4.5))
ax.plot(param.times,abs.(ξvec).^2,"g-",label="Input pulse")

ax.plot(param.times,abs.(ψ_single).^2,"ro",label="δ = 0",fillstyle="none")
#ax[2].plot(sol1.t,abs.(ref_sol)^2,"b-")
ax.plot(sol1.t,abs.(ref_sol).^2,"r-")

param.δ = 2
H = LazySum(param.δ*n,im*sqrt(param.γ/dt)*adw,-im*sqrt(param.γ/dt)*wda,param.x3/4*n*n,-param.x3/4*n)

#Define initial state
ψ_cw = onephoton(bw,ξvec)
psi = fockstate(bc,0) ⊗  ψ_cw 

#Solve
ψ = waveguide_evolution(param.times, psi, H)

#REFERENCE SOLUTION
solve_differentialeq(param,ξfun)
sol1 = solve_differentialeq(param,ξfun)
ref_sol = ξfun.(sol1.t,param.σ,param.t0)-sqrt(param.γ)*sol1

#Plot single photon waveguide state 
ψ_single = view_onephoton(ψ)/sqrt(dt)



ax.plot(param.times,abs.(ψ_single).^2,"bo",label="δ = 2",fillstyle="none")
#ax[2].plot(sol1.t,abs.(ref_sol)^2,"b-")
ax.plot(sol1.t,abs.(ref_sol).^2,"b-")


#ax[2].plot(param.times,imag.(ψ_single),"bo",fillstyle="none")

ax.set_xlabel("time [1/γ]")
#ax[2].set_xlabel("time [1/γ]")

ax.set_ylabel(L"$\xi_{out}^{(1)}$")
#ax[2].set_ylabel(L"$\xi_{out}^{(1)}$")
ax.legend()

#ax[2].set_title("Imag. part with δ = $(param.δ)")   

plt.tight_layout()