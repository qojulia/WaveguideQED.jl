using QuantumOptics
using CavityWaveguide
using LinearAlgebra
using PyPlot
pygui(true)
includet("singlecavity_simulations.jl")

#Parameters imported from singlecavity_simulations
param = parameters()
param.times = 0:0.2:20
param.t0 = 5
param.γ=2
param.δ=0

#Define operators and basis
bc = FockBasis(1)
bw = get_cwbasis(param.times,1)
a = destroy(bc)
ad = dagger(a);
n = ad*a ⊗ identityoperator(bw)
dt = param.times[2]-param.times[1]
nsteps = length(param.times)

#Precalculate hamiltonian
w = get_woper(bw,1,nsteps,1)
wda =  a ⊗ dagger(w)
adw = ad ⊗ w
H = param.δ*n + (im*sqrt(param.γ/dt))*(adw-wda)
const H_list = Array{typeof(H)}(undef,length(param.times))
for i in 1:length(param.times)
    w = get_woper(bw,1,nsteps,i)    
    wda =  a ⊗ dagger(w)
    adw = ad ⊗ w
    H_list[i] = param.δ*n + (im*sqrt(param.γ/dt))*(adw-wda)
end

function get_hamiltonian(time,psi)
    #timeindex = Int(
    return H_list[floor(Int,time/dt)+1]
end

#Define input single photon state
ξfun(t::Number,σ::Number,t0::Number) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
ξvec=sqrt(dt)*ξfun.(param.times,param.σ,param.t0)

#Define initial state
ψ_cw = Ket(bw)
ψ_cw.data[2:nsteps+1] .= ξvec
ψ0 = fockstate(bc,0) ⊗  ψ_cw 
tout, ψ = timeevolution.schroedinger_dynamic(param.times, ψ0, get_hamiltonian)

#REFERENCE SOLUTION
solve_differentialeq(param,ξfun)
sol1 = solve_differentialeq(param,ξfun)
ref_sol = ξfun.(sol1.t,param.σ,param.t0)-sqrt(param.γ)*sol1

#Plot single photon waveguide state 
ψplot = ψ[end]
ψplot = ψplot.data[1:2:end]/sqrt(dt)
ψ_single = ψplot[2:nsteps+1] 

fig,ax = subplots(1,2,figsize=(9,4.5))
ax[1].plot(param.times,real.(ψ_single))
ax[1].plot(sol1.t,real(ref_sol))
ax[2].plot(param.times,imag.(ψ_single))
ax[2].plot(sol1.t,imag(ref_sol))

ax[1].set_ylabel(L"$\xi$")
ax[1].set_title("δ = $(param.δ)")   

plt.tight_layout()
plt.savefig("plots/one_photon_contour.pdf")