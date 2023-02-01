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
#Set width of pulse
param.σ = 1

#Set simulation time
param.times = 0.0:0.1:10.0
dt = param.times[2] - param.times[1]

#Create operators for two photons interacting with cavity
bc = FockBasis(2)
bw = WaveguideBasis(2,param.times)
btotal = tensor(bc,bw)
a = sparse(destroy(bc))
ad = sparse(create(bc));
n = ad*a ⊗ identityoperator(bw)
w = destroy(bw)
wd = create(bw);
#wda = a ⊗ wd
#adw = ad ⊗ w
wda = emission(bc,bw)
adw = absorption(bc,bw)
H = param.δ*n + im*sqrt(param.γ/dt)*(adw-wda) + param.x3/4*(n*n+n)

#Define input twophoton state shape
ξfun(t1,t2,σ1,σ2,t0) = sqrt(2/σ1)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ1^2)*sqrt(2/σ2)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t2-t0)^2/σ2^2)

ψ_cw = twophoton(bw,ξfun,param.times,param.σ,param.σ,param.t0)
psi = fockstate(bc,0) ⊗ ψ_cw

ψ = waveguide_evolution(param.times, psi, H)

ψ_double = view_twophoton(ψ)
#Plot result of simulation
#Make into function to visualize?
fig,ax = subplots(1,1,figsize=(4.5,4.5))
xgrid = repeat(param.times',length(param.times),1)
ygrid = repeat(param.times,1,length(param.times))
cnt1= ax.contourf(xgrid,ygrid,ψ_double.*conj(ψ_double),100)
for c in cnt1.collections
    c.set_edgecolor("face")
end
ax.set_aspect("equal", "box")
ax.set_ylabel("time [1/γ]")
ax.set_xlabel("time [1/γ]")
ax.set_title("δ = $(param.δ),   χ3 = $(param.x3)")   
plt.tight_layout()
plt.savefig(pwd()*"/plots/two_photon_contour_xi{$(param.x3)}_delta_{$(param.δ)}.svg")