using QuantumOptics
using CavityWaveguide
using LinearAlgebra
using PyPlot
pygui(true)
includet("singlecavity_simulations.jl")
#Parameter structurer imported from singlecavity_simulations
param = parameters()
param.times = 0:0.5:10
param.t0 = 5

#Create operators for two photons interacting with cavity
bc = FockBasis(1)
bw = get_cwbasis(param.times,2)
a = destroy(bc)
ad = dagger(dense(a))
ad = sparse(ad)
n = ad*a ⊗ identityoperator(bw)
dt = param.times[2]-param.times[1]
nsteps = length(param.times)

#Precalculate operators for efficiency
w = get_woper(bw,2,nsteps,1)
wda =  a ⊗ dagger(w)
adw = ad ⊗ w
H = param.δ*n + (im*sqrt(param.γ/dt))*(adw-wda)
H_list = Array{typeof(H)}(undef,length(param.times))
for i in 1:length(param.times)
    w = get_woper(bw,2,nsteps,i)    
    wda =  a ⊗ dagger(w)
    adw = ad ⊗ w
    H_list[i] = param.δ*n + (im*sqrt(param.γ/dt))*(adw-wda)
end

function get_hamiltonian(time,psi)
    #timeindex = Int(floor(time/dt))+1
    #print("time:$time \n")
    #print("timeindex:$timeindex \n")
    return H_list[floor(Int,time/dt)+1]
end

#Define input twophoton state
ξfun(t::Number,σ::Number,t0::Number) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
ξvec=sqrt(2)*dt*tril(ξfun.(param.times,param.σ,param.t0)*transpose(ξfun.(param.times,param.σ,param.t0)))
ξvec=ξvec + (sqrt(2)/2-1)*Diagonal(ξvec)


ψ_cw = Ket(bw)
#ψ_cw.data[2:nsteps+2] .= 10
tmp = view_twophoton(ψ_cw.data,nsteps)
#tmp .= 10
tmp .= ξvec

psi = fockstate(bc,0) ⊗  ψ_cw
tout, ψ = timeevolution.schroedinger_dynamic(param.times, psi, get_hamiltonian)

ψplot = ψ[end]
#ψplot=ψ0
#timeindex = Int(floor(param.times[end]/dt)) +1
#w = get_woper(bw,2,nsteps,timeindex)    
#wda =  identityoperator(bc) ⊗ w
#ψplot = wda*wda*ψplot
ψplot = ψplot.data[1:2:end]
#ψ_single = ψplot[2:nsteps+1] 
ψ_double = view_twophoton(ψplot,nsteps)
ψ_double = ψ_double + ψ_double' - Diagonal(ψ_double)


fig,ax = subplots(1,1,figsize=(9,4.5))
xgrid = repeat(param.times',length(param.times),1)
ygrid = repeat(param.times,1,length(param.times))
ax.contourf(xgrid,ygrid,ψ_double.*conj(ψ_double),100)

ax.set_aspect("equal", "box")

ax.set_ylabel(L"$\xi$")
ax.set_title("δ = $(param.δ)")   

plt.tight_layout()
plt.savefig("plots/two_photon_contour.pdf")