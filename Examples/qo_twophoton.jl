using QuantumOptics
using CavityWaveguide
using LinearAlgebra
using PyPlot
pygui(true)
includet("singlecavity_simulations.jl")
#Parameter structure imported from singlecavity_simulations (bit overkill for now, but can prove usefull)
param = parameters()
param.δ = 0
param.γ = 1
param.t0 = 5
param.x3 = 4
param.times=0.0:0.2:10.0
dt = param.times[2] - param.times[1]
tend = param.times[end]

#Create operators for two photons interacting with cavity
bc = FockBasis(2)
bw = get_cwbasis(param.times,2)
a = sparse(destroy(bc))
ad = dagger(a);
n = ad*a ⊗ identityoperator(bw)
dt = param.times[2]-param.times[1]
nsteps = length(param.times)

#TODO put following function in module CavityWaveguide and make more general for user friendliness
#Precalculate operators for efficiency
#Can this be made to except 
function precalculate_hamiltonian(param)
    w = get_woper(bw,2,nsteps,1)
    wd = dagger(w);
    wda =  a ⊗ wd
    adw = ad ⊗ w
    H = param.δ*n + (im*sqrt(param.γ/dt))*(adw-wda) + param.x3/4*(n*n-n)
    H_list = Array{typeof(H)}(undef,length(param.times))
    for i in 1:length(param.times)
        w = get_woper(bw,2,nsteps,i)
        wda =  a ⊗ dagger(w)
        adw = ad ⊗ w
        H_list[i] = param.δ*n + (im*sqrt(param.γ/dt))*(adw-wda)+ param.x3/4*(n*n-n)
    end
    return H_list
end

function get_hamiltonian(time,psi)
    return H_list[floor(Int,time/dt)+1]
end

function fout(time,psi)
    if time == tend
        return psi
    else
        return 0
    end
end

#Define input twophoton state shape
#Can this be done in a nicer way?
ξfun(t::Number,σ::Number,t0::Number) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
ξvec=sqrt(2)*dt*tril(ξfun.(param.times,param.σ,param.t0)*transpose(ξfun.(param.times,param.σ,param.t0)))

#Define input waveguide state
ψ_cw = Ket(bw)
tmp = view_twophoton(ψ_cw.data,nsteps)
tmp .= ξvec

#Tensor with cavity
psi = fockstate(bc,0) ⊗ ψ_cw

#Precalc hamiltonian
H_list = precalculate_hamiltonian(param)

#Solve
tout, ψ = timeevolution.schroedinger_dynamic(param.times, psi, get_hamiltonian,fout=fout)

#Extract last vector
ψplot = ψ[end].data

#TODO: Make the following into function that views specified output
ψplot = ψplot[1:3:end]
ψ_double = view_twophoton(ψplot,nsteps)
ψ_single = ψplot[2:nsteps+1] 
ψ_double = ψ_double + ψ_double' - Diagonal(ψ_double)

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