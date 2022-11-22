using CavityWaveguide
using Interpolations
using PyPlot
pygui(true)

includet("singlecavity_simulations.jl")

param = parameters()
param.δ = 1
param.times=0.0:0.05:50.0
@unpack times = param
dt = times[2] - times[1]

ξfun(t::Number,σ::Number,t0::Number) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
ξvec = ξfun.(param.times,param.σ,param.t0)

ξout = solve_onephotonwaveguide(param,ξvec)
sol1 = solve_differentialeq(param,ξfun)
ref_sol = ξfun.(sol1.t,param.σ,param.t0)-sqrt(param.γ)*sol1

fig,ax = subplots(1,1,figsize=(9,7))
ax.plot(times,real.(ξout.data[2].ξ01),"ro",label="Time bin",fillstyle="none")
ax.plot(sol1.t,imag.(ref_sol),"b-",label="Diff. Eq.")

ax.plot(sol1.t,real.(ref_sol),"r-",label="Real part")
ax.plot(times,imag.(ξout.data[2].ξ01),"bo",fillstyle="none",label="Imag. part")



ax.set_ylabel(L"$\xi$")
ax.set_title("δ = $(param.δ)")

ax.legend(loc="best")
#ax[2].set_xlabel("time")

interp_new = LinearInterpolation(times,ξout.data[2].ξ01)

#ax[2].semilogy(sol1.t,abs.(real.(ref_sol - interp_new(sol1.t))./(real.(interp_new(sol1.t)))),"ro")
#ax[2].semilogy(sol1.t,abs.(imag.(ref_sol - interp_new(sol1.t))./(imag.(interp_new(sol1.t)))),"bo")

#ax[2].semilogy(sol1.t,abs.(ξ.(sol1.t,param1.σ,param1.t0)-sqrt(param1.γ)*sol1-interp_old(sol1.t))./(ξ.(sol1.t,param1.σ,param1.t0)-sqrt(param1.γ)*sol1),"ro")
#ax[2].set_ylim(1e-5,0.01)
plt.tight_layout()
plt.savefig("plots/diff_vs_timebin_v2.pdf")