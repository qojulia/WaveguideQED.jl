using QuantumOptics
using CavityWaveguide
using LinearAlgebra
using PyPlot
pygui(true)
include("singlecavity_simulations.jl")
#Parameter structure imported from singlecavity_simulations (bit overkill for now, but can prove usefull)
param = parameters()
param.δ = 0
param.γ = 1
param.t0 = 5
#param.x3 = 4
param.times=0.0:0.2:20.0
dt = param.times[2] - param.times[1]
tend = param.times[end]

#Create operators for two photons interacting with cavity
bc = FockBasis(1)
bw = WaveguideBasis(1,param.times)
btotal = tensor(bc,bw)
a = sparse(destroy(bc))
ad = sparse(create(bc));
n = ad*a ⊗ identityoperator(bw)
dt = param.times[2]-param.times[1]

nsteps = length(param.times)

n = ad*a ⊗ identityoperator(bw)
dt = param.times[2]-param.times[1]
nsteps = length(param.times)
w = destroy(bw)
wd = create(bw);
wda = LazyTensor(btotal,btotal,[1,2],(a,wd))
adw = LazyTensor(btotal,btotal,[1,2],(ad,w))
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

fig,ax = subplots(1,2,figsize=(9,4.5))
ax[1].plot(param.times,real.(ψ_single))
ax[1].plot(sol1.t,real(ref_sol))
ax[2].plot(param.times,imag.(ψ_single))
ax[2].plot(sol1.t,imag(ref_sol))

ax[1].set_ylabel(L"$\xi$")
ax[1].set_title("δ = $(param.δ)")   

plt.tight_layout()
plt.savefig("plots/one_photon_contour.pdf")