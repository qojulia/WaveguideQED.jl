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
param.x3 = 4
param.times=0.0:0.1:10.0
dt = param.times[2] - param.times[1]
tend = param.times[end]

#Create operators for two photons interacting with cavity
bc = FockBasis(2)
bw = WaveguideBasis(2,param.times)
btotal = tensor(bc,bw)
a = sparse(destroy(bc))
ad = sparse(create(bc));

n = ad*a ⊗ identityoperator(bw)
dt = param.times[2]-param.times[1]
nsteps = length(param.times)
w = destroy(bw)
wd = create(bw);
wda = LazyTensor(btotal,btotal,[1,2],(a,wd))
adw = LazyTensor(btotal,btotal,[1,2],(ad,w))
H = LazySum(param.δ*n,im*sqrt(param.γ/dt)*adw,-im*sqrt(param.γ/dt)*wda,param.x3/4*n*n,-param.x3/4*n)

#Define input twophoton state shape
#Can this be done in a nicer way?
ξfun(t::Number,σ::Number,t0::Number) = 1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2)/sqrt(0.2820947917738782)
ξvec=sqrt(2)*dt*tril(ξfun.(param.times,param.σ,param.t0)*transpose(ξfun.(param.times,param.σ,param.t0)))

#Define input waveguide state
ψ_cw = Ket(bw)
tmp = view_twophoton(ψ_cw)
tmp .= ξvec

#Tensor with cavity
psi = fockstate(bc,0) ⊗ ψ_cw

#Solve
ψ = waveguide_evolution(param.times, psi, H)

#TODO: Make the following into function that views specified output
ψ_double = view_twophoton(ψ)
ψ_double = ψ_double + ψ_double' - Diagonal(ψ_double)
ψ_single = view_singlephoton(ψ)

#Plot result of simulation
#Make into function to visualize?
fig,ax = subplots(1,1,figsize=(9,4.5))
xgrid = repeat(param.times',length(param.times),1)
ygrid = repeat(param.times,1,length(param.times))
ax.contourf(xgrid,ygrid,ψ_double.*conj(ψ_double),100)
ax.set_aspect("equal", "box")
ax.set_ylabel(L"$\xi$")
ax.set_title("δ = $(param.δ)")   
plt.tight_layout()
plt.savefig("plots/two_photon_contour.pdf")