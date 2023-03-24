# Input-Output Waveguides

In many of the examples considered so far, we only consider a single waveguide that serves as both input and output, thus only allowing for only one-sided cavities or quantum systems at the end of a waveguide. A more realistic scenario is having a waveguide with a quantum system in the middle. Here an incoming waveguide carying an excitation could scatter on the quantum system and one would have excitations going away from the quantum system in both the first and latter part of the waveguide as illustrated here:

![`alt text`](two_waveguide_lodahl.png)

A way to model this scenario is to have to waveguides: a waveguide to the left and the right, describing the first half of the waveguide and the latter part of the waveguide. For this we use the custom basis called [`LeftRightWaveguideBasis`](@ref) (see [`Two Waveguides`](@ref twowaveguide) for an introduction). We initialize [`LeftRightWaveguideBasis`](@ref) and a basis for the atom:

```jldoctest
bw = LeftRightWaveguideBasis(2,times)
be = FockBasis(1)
```

When creating operators, we now have to specify which waveguide they are acting on (input or output). The equivalent of [`emission`](@ref) and [`absorption`](@ref) is then:

```jldoctest
wdL = leftemission(bw,be)
wL = leftabsorption(bw,be)
wdR = rightemission(bw,be) 
wR = rightabsorption(bw,be)
```

where $wdL = $`


```jldoctest
using CavityWaveguide
using QuantumOptics
using PyPlot
pygui(false)
pygui(true)

times = 0:0.1:10
dt = times[2] - times[1]
κ1 = 1
κ2 = 1

#Create operators for two    photons interacting with cavity
be = FockBasis(1)
bw = LeftRightWaveguideBasis(2,times)
bw_single = LeftRightWaveguideBasis(1,times)
b = bw ⊗ be
b_single = bw_single ⊗ be
a = embed(b,2,destroy(be))
ad = embed(b,2,create(be))
wdL = leftemission(bw,be)
wL = leftabsorption(bw,be)
wdR = rightemission(bw,be) 
wR = rightabsorption(bw,be)

a_single = embed(b_single,2,destroy(be))
ad_single = embed(b_single,2,create(be))
wdL_single = leftemission(bw_single,be)
wL_single = leftabsorption(bw_single,be)
wdR_single = rightemission(bw_single,be) 
wR_single = rightabsorption(bw_single,be)

H = im*sqrt(κ1/dt)*(wL-wdL) + im*sqrt(κ2/dt)*(wR-wdR)
H_single = im*sqrt(κ1/dt)*(wL_single-wdL_single) + im*sqrt(κ2/dt)*(wR_single-wdR_single) 


ξfun(t1,t2,σ1,σ2,t0) = sqrt(2/σ1)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ1^2)*sqrt(2/σ2)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t2-t0)^2/σ2^2)
ξ_one_fun(t1,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ^2)

psi_double_list = []
psi_single_list = []
widths = [0.8,1,1.5,2.9,8.3]/sqrt(4*log(2))
t0 = 5

for w in widths
    psi_in = twophoton(bw,:input,ξfun,times,w,w,t0) ⊗ fockstate(be,0)
    psi_in_single = onephoton(bw_single,:input,ξ_one_fun,times,w,t0) ⊗ fockstate(be,0)
    psi_out = waveguide_evolution(times,psi_in,H)
    psi_out_single = waveguide_evolution(times,psi_in_single,H_single)
    psi_R_scat = TwoPhotonView(psi_out,type=:input)
    push!(psi_double_list,psi_R_scat)

    psi_R_scat_single = zeros(ComplexF64,(length(times),length(times)))
    single_R = OnePhotonView(psi_out_single,type=:input)

    for i in eachindex(times)
        for j in eachindex(times)
            psi_R_scat_single[i,j] = single_R[i]*single_R[j]
        end
    end
    push!(psi_single_list,psi_R_scat_single)
end

fig,axs = subplots(2,5,figsize=(18,7))
times = times .- 5
xgrid = repeat(times',length(times),1)
ygrid = repeat(times,1,length(times))
for (i,ax) in enumerate(axs[2,:])
    cnt1= ax.contourf(xgrid,ygrid,psi_double_list[i].*conj(psi_double_list[i]),100,cmap="Blues")
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    ax.set_xlabel(L"$t_1$")
    ax.set_ylabel(L"$t_2$")
    ax.set_title("")
end

for (i,ax) in enumerate(axs[1,:])
    cnt1= ax.contourf(xgrid,ygrid,psi_single_list[i].*conj(psi_single_list[i]),100,cmap="Reds")
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    ax.set_yticks([])
    ax.set_xlabel(L"$t_1$")
end
plt.tight_layout()
plt.savefig(pwd()*"/plots/lodahl_fig2.pdf")
```