using QuantumOptics
using CavityWaveguide
using LinearAlgebra
using PyPlot
pygui(true)
include("singlecavity_simulations.jl")
#Parameter structure imported from singlecavity_simulations (bit overkill for now, but can prove usefull)
param = parameters()

#Set detuning:
param.δ = 0
#Set non linearity:
param.x3 = 0

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
H = LazySum(param.δ*n,im*sqrt(param.γ/dt)*adw,-im*sqrt(param.γ/dt)*wda,param.x3/4*n*n,-param.x3/4*n)


#Define input twophoton state shape
ξfun(t::Number,σ::Number,t0::Number) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
ξvec=sqrt(dt)*ξfun.(param.times,param.σ,param.t0)

#Define initial state
ψ_cw = Ket(bw)
tmp = view_singlephoton(ψ_cw)
tmp .= ξvec
psi = fockstate(bc,0) ⊗  ψ_cw 

#Solve
ψ = waveguide_evolution(param.times, psi, H)

#REFERENCE SOLUTION
solve_differentialeq(param,ξfun)
sol1 = solve_differentialeq(param,ξfun)
ref_sol = ξfun.(sol1.t,param.σ,param.t0)-sqrt(param.γ)*sol1

#Plot single photon waveguide state 
ψ_single = view_singlephoton(ψ)/sqrt(dt)

fig,ax = subplots(2,1,figsize=(9,9))
ax[1].plot(param.times,real.(ψ_single),"ro",label="CavityWaveguide.jl",fillstyle="none")
ax[2].plot(sol1.t,imag.(ref_sol),"b-")
ax[1].plot(sol1.t,real.(ref_sol),"r-",label="Reference sol.")
ax[2].plot(param.times,imag.(ψ_single),"bo",fillstyle="none")

ax[1].set_xlabel("time [1/γ]")
ax[2].set_xlabel("time [1/γ]")

ax[1].set_ylabel(L"$\xi_{out}^{(1)}$")
ax[2].set_ylabel(L"$\xi_{out}^{(1)}$")
ax[1].legend()

ax[1].set_title("Real part with δ = $(param.δ)")   
ax[2].set_title("Imag. part with δ = $(param.δ)")   

plt.tight_layout()
