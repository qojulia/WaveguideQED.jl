using CavityWaveguide
using Interpolations
using PyPlot
using LinearAlgebra
pygui(true)

includet("singlecavity_simulations.jl")

param = parameters()
param.δ = 0
param.x3 = 0
param.γ = 1
param.t0 = 5
param.times=0.0:0.05:20.0
@unpack times = param
dt = times[2] - times[1]

ξfun(t::Number,σ::Number,t0::Number) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
ξvec=sqrt(2)*dt*tril(ξfun.(times,param.σ,param.t0)*transpose(ξfun.(times,param.σ,param.t0)))
ξvec=ξvec + (sqrt(2)/2-1)*Diagonal(ξvec)


ξout = solve_twophotonwaveguide_kerr(param,ξvec)

a = CavityWaveguide.annihilation(ξout.basis)
ad = CavityWaveguide.creation(ξout)
n  = CavityWaveguide.number(ξout)

fig,ax = subplots(1,2,figsize=(9,4.5))
xgrid = repeat(times',length(times),1)
ygrid = repeat(times,1,length(times))
cnt1 =  ax[1].contourf(xgrid,ygrid,ξout.data[3].ξ02 .* conj.(ξout.data[3].ξ02),100)
ax[1].set_aspect("equal", "box")
#ax.contourf(xgrid,ygrid,real(ξ2))
#Z = [x*y for x in reverse(times),y in times]
#
#ax[i].plot(time[i],ξout.ξ1./output.-1)
#ax[1].set_ylabel(L"$|\xi_{out}^{(2)}|^2$")
ax[1].set_title(L"$\chi_3$ = "*"$(param.x3)")
ax[1].set_xlabel(L"t [$1/\gamma$]")
ax[1].set_ylabel(L"t [$1/\gamma$]")


param.δ = 0
param.x3 = 4
ξvec=tril(ξfun.(times,param.σ,param.t0)*transpose(ξfun.(times,param.σ,param.t0))*sqrt(2)*dt)
ξvec=ξvec + (sqrt(2)/2-1)*Diagonal(ξvec)

ξout = solve_twophotonwaveguide_kerr(param,ξvec)

cnt2 = ax[2].contourf(xgrid,ygrid,ξout.data[3].ξ02 .* conj.(ξout.data[3].ξ02),100)

ax[2].set_aspect("equal", "box")
ax[2].set_xlabel(L"t [$1/\gamma$]")
ax[2].set_ylabel(L"t [$1/\gamma$]")

ax[2].set_title(L"$\chi_3$ = "*"$(param.x3)")

for c in cnt1.collections
    c.set_edgecolor("face")
end
for c in cnt2.collections
    c.set_edgecolor("face")
end

#ax.legend(loc="best")
#ax[2].plot(times,real.(ξout.data[3].ξ11))
#ax[2].plot(times,imag.(ξout.data[3].ξ11))
plt.tight_layout()
plt.savefig("plots/two_photon_contour_kerr.pdf")
