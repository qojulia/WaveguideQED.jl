# [Two Waveguides](@id twowaveguide)
In the previous examples, we have only considered cases with a single waveguide. In this toturial, we show how to model a beamsplitter and an optical switch using two waveguides. A beamsplitter or a swap gate can be modelled using the Hamiltonian $H = V(w_a^\dagger w_b + w_b^\dagger w_a)$ where V is some interaction strength that determines which interaction is moddeled (we will discuss this in detail later). $w_a$ and $w_b$ is the annihilation operators of the two waveguides. We can describe the state of two waveguides with a total of N excitations by adding an argument specifying the number of waveguides as:

```jldoctest
times = 0:0.1:10
NPhotons = 2
NWaveguides = 2
bw_twophotons = WaveguideBasis(NPhotons,NWaveguides,times)
```

When creating operators, we now have to specify which waveguide they are acting on (in this case number one or two). This is done by an extra argument to [`create`](@ref) and [`destroy`](@ref):

```jldoctest
wdL = create(bw,1)
wL = destroy(bw,1)
wdR = create(bw,2) 
wR = destroy(bw,2)
```

Similarly, initializing one or two photon states in the first or second waveguide is done by:

```jldoctest
ξ(t,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2)
ξ2(t1,t2,σ,t0) = ξ(t1,σ,t0)*ξ(t2,σ,t0)
ψ_single_first = onephoton(bw,1,ξ,times,2,5)
ψ_double_first = twophoton(bw,1,ξ2,times,2,5)
ψ_single_second = onephoton(bw,2,ξ,times,2,5)
ψ_double_second = twophoton(bw,2,ξ2,times,2,5)
```

If we want to describe a simultanous excitation in both waveguides (states like $\ket{1_i}_\mathrm{left}\ket{1_j }_\mathrm{right}$) we specify both indeces of the waveguides:

```jldoctest
ψ_single_first_and_second = twophoton(bw,[1,2],ξ2,times,2,5)
```

## Beamsplitter
Let's now treat the same example as in [Interference on Beamsplitter](@ref BStoturial). We consider the two waveguides in a identic single photon state and thus use the above defined `ψ_single_first_and_second`. The Hamiltonian governing a beamsplitter in the time binned formalism has $V= \pi/4$:

```jldoctest
V = pi/4
H = im*V/dt*(wdR*wL - wdL*wR)
```

We can then evolve the system under this Hamiltonian to perform the beamsplitting operation:

```jldoctest
psi_out = waveguide_evolution(times,ψ_single_first_and_second,H)
```

We can then view the final state to verify that we only have twophotons in the same waveguide simultanouesly:

```jldoctest
psi_second = TwoPhotonView(psi_out,2)
psi_first = TwoPhotonView(psi_out,1)
psi_first_second = TwoPhotonView(psi_out,[1,2])
julia> norm(psi_R)^2
julia> norm(psi_L)^2
julia> norm(psi_LR)^2
0.49999981822067935
0.49999981822067935
8.736388404016349e-7
```

Except for numerical errors we thus have 50% chance of observing both photons in the same waveguide and 0 (8.736388404016349e-9)% of observing both photons in each of the waveguide simultanoues. 

## Swap
If we instead choose $V = \pi / 2$ we get the SWAP operation. Let us consider on photons in the left waveguide and swap them to right waveguide and plot before and after:

```jldoctest
V = pi/2
H = im*V/dt*(wdR*wL - wdL*wR)
psi_out_swap = waveguide_evolution(times,ψ_single_first,H)
first_before = OnePhotonView(ψ_single_first,1)
second_before = OnePhotonView(ψ_single_first,2)
first_after = OnePhotonView(psi_out_swap,1)
second_after = OnePhotonView(psi_out_swap,2)

fig,ax = subplots(1,2,figsize=(9,4.5))
ax[1].plot(times,first_before,"r-",label="First")
ax[1].plot(times,second_before,"b-",label="Second")
ax[2].plot(times,first_after,"r-")
ax[2].plot(times,second_after,"b-")

ax[1].legend(loc="lower left")
ax[1].set_title("Before")
ax[1].set_xlabel("time [a.u]")
ax[1].set_ylabel("ξ(t)")

ax[2].set_title("After")
ax[2].set_xlabel("time [a.u]")

plt.tight_layout()
```
!["Alt text"](swap.png)



!!! info "WaveguideBasis(2,2,times) vs. Waveguide $\otimes$ Waveguide"
    Instead of using the custom basis for handling two waveguides, one could instead just do a tensor product between two waveguides basises. This naive approach would look something like:
    
    ```jldoctest
    times = 0:0.1:10
    bw = WaveguideBasis(2,times)
    Btotal = bw ⊗ bw
    ```

    This might work if you only consider single photon excitations in the waveguides, but if you go consider two photon excitations (as in the above) the hilbert space blows up. Indeed, in the above example, the hilbert space is of size: 27.594.009!!! However, since we often know that the system in total only has two excitations, there is no possibility of having two photons in both waveguides simultanously (states of type: $\ket{1_k,1_j}_{left}\ket{1_l,1_m}_{right}$). This part of the Hilber space takes up the majority since it scales as $\propto N^4$ where N is the number of timebins. Instead we can exploit that only a total of two excitations is present simultanosly in the system. For this we use the custom basis [`WaveguideBasis`](@ref) with an addination input stating the number of waveguides:

    ```jldoctest
    bw = WaveguideBasis(2,2,times)
    ```
    The hilbert space is now of size: 20706, this is of course still large Hilbert space, but 3 orders of magnitude smaller than the naive approach as it still only scales as $\propto N^2$. 

