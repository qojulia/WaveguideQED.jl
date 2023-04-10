# Scattering on two level system
In the following, we show that with our framwework we can reproduce the theoretical results obtained in [Le Jeannic, et al. Nat. Phys. 18, 1191–1195 (2022)](https://www.nature.com/articles/s41567-022-01720-x) 


In many of the examples considered so far, we only consider a single waveguide that serves as both input and output, thus only allowing for only one-sided cavities or quantum systems at the end of a waveguide. A more realistic scenario is having a waveguide with a quantum system in the middle. Here an incoming waveguide carying an excitation could scatter on the quantum system and one would have excitations going away from the quantum system in both the first and latter part of the waveguide as illustrated here:[^1]

[^1]: [DynamicalPhotonLodahl2022](@cite)

![`alt text`](./illustrations/two_waveguide_lodahl.png)


A way to model this scenario is to have two waveguides: a waveguide to the left and the right, describing the first half of the waveguide and the latter part of the waveguide. For this we use [`WaveguideBasis`](@ref) but with an extra argument specifying that we need 2 waveguides (see [`Two Waveguides`](@ref twowaveguide) for an introduction). We initialize [`WaveguideBasis`](@ref) with two waveguides and a basis for the atom (note that a fockbasis with only one excitation allowed is the same as a two-level-system):

```@example lodahl
using WaveguideQED #hide
using QuantumOptics #hide
times = 0:0.1:10
dt = times[2] - times[1]
bw = WaveguideBasis(2,2,times)
be = FockBasis(1)
nothing #hide
```

We then define the operators for the interaction between atom and waveguide as (notice the second argument in create(bw,1) that defines which waveguide we are adressing):

```@example lodahl
wdLa = create(bw,1) ⊗ destroy(be)
adwL = destroy(bw,1) ⊗ create(be)
wdRa = create(bw,2) ⊗ destroy(be)
adwR = destroy(bw,2) ⊗ create(be)
nothing #hide
```

where $\mathrm{wdLa} = w_L ^\dagger a$, $\mathrm{wdRa} = w_R ^\dagger a$, $\mathrm{adwL} = w_L  a^\dagger$, and $\mathrm{adwR} = w_R  a^\dagger$. In this example, we, however, also need an interaction between the waveguides. We therefore we define the creation and annihilation operators 

```@example lodahl
wdL = create(bw,1) ⊗ identityoperator(be)
wL = destroy(bw,1) ⊗ identityoperator(be)
wdR = create(bw,2) ⊗ identityoperator(be)
wR = destroy(bw,2) ⊗ identityoperator(be)
nothing #hide
```

The interaction should carry over the momentum of the left pulse into the right waveguide and the interaction should therefore model a SWAP gate. This corresponds to $V = \pi /2$ and thus we have the interaction Hamiltonian:

```@example lodahl
V = pi/2
κ1 = 1
κ2 = 1
H = im*sqrt(κ1/dt)*(adwL-wdLa) + im*sqrt(κ2/dt)*(wdRa-adwR) + V/dt *(wdR*wL + wdL* wR)
nothing #hide
```

We can now study how single or two photon states scatter on the atom. We define the initial onephoton or twophoton gaussian state and solve it using the defined Hamiltonian:


```@example lodahl
ξ₁(t1,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ^2)
ξ₂(t1,t2,σ1,σ2,t0) = ξ₁(t1,σ1,t0) * ξ₁(t2,σ2,t0) 
w = 1
t0 = 5
ψ1 = onephoton(bw,1,ξ₁,times,w,t0) ⊗ fockstate(be,0)
ψ2 = twophoton(bw,1,ξ₂,times,w,w,t0) ⊗ fockstate(be,0)
ψScat1 = waveguide_evolution(times,ψ1,H)
ψScat2 = waveguide_evolution(times,ψ2,H)
nothing #hide
```

Viewing the scattered states is then done using TwoPhotonView and the index for the corresponding waveguide. Giving two indeces returns instead the combined single photon state in both waveguides $\sum_{j,k} \ket{1_j}_1 \ket{1_k}_2$:

```@example lodahl
ψ2LeftScat = TwoPhotonView(ψScat2,[:,1],1)
ψ2RightScat = TwoPhotonView(ψScat2,[:,1],2)
ψ2LeftRightScat = TwoPhotonView(ψScat2,[:,1],2,1)
nothing #hide
```

For the single photon states we have to calculate the two time scattered distribution as:

```@example lodahl
ψ1LeftScat = zeros(ComplexF64,(length(times),length(times)))
ψ1RightScat = zeros(ComplexF64,(length(times),length(times)))
ψ1LeftRightScat = zeros(ComplexF64,(length(times),length(times)))
ψ1Right = OnePhotonView(ψScat1,[:,1],1)
ψ1Left = OnePhotonView(ψScat1,[:,1],2)

for i in eachindex(times)
    for j in eachindex(times)
        ψ1LeftScat[i,j] = ψ1Left[i]*ψ1Left[j]
        ψ1RightScat[i,j] = ψ1Right[i]*ψ1Right[j]
        ψ1LeftRightScat[i,j] = ψ1Left[i]*ψ1Right[j]
    end
end
nothing #hide
```

Finall, this can be plotted and we note that this matches fig. 3 in Ref. ^[1]:

```@example lodahl
using PyPlot; #hide
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams"); #hide
rcParams["font.size"] = 20; #hide
fig,axs = subplots(3,2,figsize=(9,17))
plot_list = [ψ2LeftScat,ψ2RightScat,ψ2LeftRightScat,ψ1LeftScat,ψ1RightScat,ψ1LeftRightScat]
for (i,ax) in enumerate(axs)
    plot_twophoton!(ax,plot_list[i],times)
end
axs[1].set_ylabel("\$C^{RR}\$ \n t2 [a.u]")
axs[2].set_ylabel("\$C^{LL}\$ \n t2 [a.u]")
axs[3].set_ylabel("\$C^{LR}\$ \n t2 [a.u]")
axs[3].set_xlabel("t1 [a.u]")
axs[6].set_xlabel("t1 [a.u]")
plt.tight_layout()
plt.savefig("lodahl_fig3.svg") #hide
nothing #hide
```
![lodahl](lodahl_fig3.svg)

If we consider the single photon state, we can also visualize the temporal evolution as:

![alt text](./animations/lodahl_onephoton_gif.gif)


