# Theoretical Background
In this section, we go over the necessary theory to work with continous fockstates in the **WaveguideQED.jl**


## Continous Fock States

The single photon continuous fock state can be defined as:

$$\begin{equation*}
    \ket{\psi} = W^\dagger(\xi) \ket{0} = \int_{t_0}^{t_{end}} \mathrm{d}t \ \xi(t) w^\dagger(t) \ket{\emptyset}
\end{equation*}$$

here $W^\dagger(\xi)$ creates a photon with the wavefunction $\xi(t)$. $w^\dagger(t)$ is the creation operator for a photon at time $t$, and it obeys the commutation relation: $\left[w(t),w(t')\right ] = \delta(t-t')$. The probability of observing a photon at time $t$ is given by: $\bra{0} w(t) \ket{\psi} = |\xi(t)|^2$. The wavefunction $\xi(t)$ thus describes the temporal distribution of the photon.

The heart of the photon time binning is discretizing the continuous fock state into time bins of width $\Delta t$. The interaction with the emitter/cavity will is then assumed to span only one timebin at a time, corresponding to a spectrally flat interaction between the waveguide and emitter/cavity. We thus discretize the annihilation and creation operators by taking[^1]:

$$\begin{equation*}
    w(t_k) = w(k \Delta t) \rightarrow  \frac{w_k}{\sqrt{\Delta t}} \ \ \  \text{with} \ \left[ w_j, w_k^\dagger \right ] = \delta_{jk}
\end{equation*}$$

where $w_k$ is the descritized operator and the factor of $1/\sqrt{\Delta t}$ assures the commutator relation in the limit of $\Delta t \rightarrow 0$. We denote the action of the discretized creation operator as: $w_k\dagger \ket{\emptyset} = \ket{1_k}$ meaning a single photon in timebin $k$. This means that the single photon continuous fock state becomes:

$$\begin{equation*}
    \ket{\psi} = \int_{t_0}^{t_{end}} \mathrm{d}t \ \xi(t) w^\dagger(t) \ket{\emptyset} \rightarrow 
\sum_{k=1}^N \sqrt{\Delta t} \xi(t_k) w_k^\dagger \ket{\emptyset}
\end{equation*}$$

In `WaveguideQED.jl`, the timebins above are represented as elements in arrays corresponding to each timebin:

![Alt text](./illustrations/onephoton_array.png)

Let`s say you want to represent a single photon contionous fock state that starts at $t=0$ and ends at $t=10$ with $\Delta t = 0.1$. This can be done be first creating waveguide basis defined on such a timeinterval:

```@example theory
using WaveguideQED
times = 0:0.1:10
bw = WaveguideBasis(1,times)
nothing #hide
```

Notice that the input for WaveguideBasis is `1` and `times`. `1` denotes the maximum excitation number of fockstates (currently can only be 1 or 2) and `times` the time interval over which the continous fock state is defined. To define the continous fockstate we need to give a wavefuntion $\xi$. In the following we define a gaussian wavefunction located around $t=5$ with a width of $\sigma = 1$:

```@example theory
ξ(t,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2)
σ,t0 = 1,5
ψ = onephoton(bw,ξ,times,σ,t0)
nothing #hide
```

This state can be visuallized by:

```@example theory
using PyPlot
viewed_state = OnePhotonView(ψ)
fig,ax = subplots(1,1,figsize=(9,4.5))
ax.plot(times,viewed_state,"r-")
ax.set_xlabel("Time [a.u]")
ax.set_ylabel(L"$\xi(t)$")
plt.tight_layout()
plt.savefig("one_continous_fockstate.svg") #hide
nothing #hide
```
![alt text](one_continous_fockstate.svg)

The timebinned creation and annihilation operators are easily created from the basis:

```@example theory
w = destroy(bw)
wd = create(bw)
nothing #hide
```

The timebin that the operator acts on is set by either:

```@example theory
w.timeindex = 10
wd.timeindex = 10
nothing #hide
```
 or:

```@example theory
set_waveguidetimeindex!(w,10)
set_waveguidetimeindex!(wd,10)
nothing #hide
```

The effect of the creation operator is to create a photon in timebin k and can be illustrated as:

![Alt text](./illustrations/one_photon_creation.png)

This is also seen if we plot the creation operator acting on the vacuum:

```@example theory
using PyPlot
ψ = wd*zerophoton(bw)
viewed_state = OnePhotonView(ψ)
fig,ax = subplots(1,1,figsize=(9,4.5))
ax.plot(times,viewed_state,"r-")
ax.set_xlabel("Time [a.u]")
ax.set_ylabel(L"$\xi(t)$")
plt.tight_layout()
plt.savefig("created_onephoton_continous_fockstate.svg") #hide
nothing #hide
```
![alt text](created_onephoton_continous_fockstate.svg)

We see a spike around `t = times[10] = 0.9` where we now created an excitation. In itself, the waveguide basis, states, and operators are not particularly interesting, but when combined with other quantum mechanical systems such as cavities and emitters the framework can produce powerfull results. See the next example of how to combine with quantum systems defined in ['QuantumOptics.jl'](https://qojulia.org/)


## Continous two-photon fock states

So far we have considered only one excitation in the waveguide. We can extend the definition of a one-photon continous fock state to a two-photon state as[^2]:

$$\begin{align*}
\frac{1}{\sqrt{2}}\left[W^\dagger(\xi)\right]^2|0\rangle &= \frac{1}{\sqrt{2}} \int_{t_0}^{t_{end}} d t^{\prime} \int_{t_0}^{t_{end}} d t \ \xi(t) \xi\left(t^{\prime}\right) w^\dagger(t) w^\dagger\left(t^{\prime}\right)|0\rangle \\
 &= \frac{1}{\sqrt{2}} \int_{t_0}^{t_{end}} d t^{\prime} \int_{t_0}^{t_{end}} d t \ \xi^{(2)}(t,t') w^\dagger(t) w^\dagger\left(t^{\prime}\right)|0\rangle  
\end{align*}$$

Here, we here defined the two photon wavefunction $$\xi^{(2)}(t,t') = \xi(t) \xi\left(t^{\prime}\right)$$. The state is now defined over two times, which describes the probability of observing photon A at time $t$ and photon B at time $t'$. In this case, the state is a product state $$\xi^{(2)}(t,t') = \xi(t) \xi\left(t^{\prime}\right)$$, and both probabilities are described by the (same) single photon wavefunction $\xi(t)$, but one could have entangled states across time. This means a non-seperable wavefunction $$\xi^{(2)}(t,t') \neq \xi_1(t)\xi_2(t')$$. For now we will consider a symmetric and seperable state.

The timebinning is in a similar fashion defined as:

$$\begin{align*}
\frac{1}{\sqrt{2}}\left[W^\dagger(\xi)\right]^2|0\rangle &= \frac{1}{\sqrt{2}} \int_{t_0}^{t_{end}} d t^{\prime} \int_{t_0}^{t_{end}} d t \ \xi(t) \xi\left(t^{\prime}\right) w^\dagger(t) w^\dagger\left(t^{\prime}\right)|0\rangle \\
& \rightarrow \frac{1}{\sqrt{2}} \sum_{i=1}^N \sum_{k=1}^N \xi\left(t_i\right) \xi\left(t_k\right) w^\dagger\left(t_i\right) w^{\dagger}\left(t_k\right)|0\rangle \\
& =\frac{1}{\sqrt{2}} \sum_{i=1}^N \sum_{k \neq i}^N \xi\left(t_i\right) \xi\left(t_k\right) w^{\dagger}\left(t_i\right) w^{\dagger}\left(t_k\right)|0\rangle+\sum_{i=1}^N \xi\left(t_i\right) \xi\left(t_i\right)\left|2 t_i\right\rangle \\
& =\frac{2}{\sqrt{2}} \sum_{i=1}^N \sum_{k>i}^N \xi\left(t_i\right) \xi\left(t_k\right)\left|1_{t_i} 1_{t_k}\right\rangle+\sum_{i=1}^N \xi\left(t_i\right) \xi\left(t_i\right)\left|2 t_i\right\rangle \\
& =\sqrt{2} \sum_{i=1}^N \sum_{k > i}^N \xi\left(t_i\right) \xi\left(t_k\right) \mid 1_{t_i} 1_{t_k}\rangle + \sum_{i=1}^N \xi\left(t_i\right) \xi\left(t_i\right)\left|2 t_i\right\rangle
\end{align*}$$


The sum is allowed to run over only half of the times due to the symmetry of the photons (it's equivalent having one photon at time bin k and then one photon at time bin j or one photon at time bin j and then one photon at time bin k). This is how the twophoton state is saved in the underlying arrays and can be illustrated as:

![alt text](./illustrations/twophoton_array.png)

Creating then is:

![alt text](./illustrations/two_photon_creation.png)

We can define a twophoton basis and corresponding operator by:

```@example theory
bw = WaveguideBasis(2,times)
wd = create(bw)
```
The creation operator can then be visualized by acting on a onephoton state with ones all over. This is seen in the following. Note that the state is visualized as a contour plot mirrored around the diagonal.

```@example theory
set_waveguidetimeindex!(wd,50)
psi_plot = wd*onephoton(bw,x->1,times)
fig,ax = subplots(1,1,figsize=(4.5,4.5))
plot_twophoton!(ax,psi_plot,times)
plt.savefig("twophoton_created.svg") #hide
nothing #hide
```

![alt text](twophoton_created.svg)

If we want to create a twophotn gaussian state we instead do:

```@example theory
ξ(t,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2)
ξ2(t1,t2,σ,t0) = ξ(t1,σ,t0)*ξ(t2,σ,t0)
σ,t0 = 1,5
ψ = twophoton(bw,ξ2,times,σ,t0)
nothing #hide
```

Here, we defined the two photon equivalent of our single photon gaussian state. When we visualize it, we now need two times, and we make a contour plot. This is easily done viewing the twophoton state: 

```@example theory
using PyPlot
viewed_state = TwoPhotonView(ψ)
fig,ax = subplots(1,1,figsize=(4.5,4.5))
plot_twophoton!(ax,viewed_state,times)
ax.set_ylabel("time [1/γ]")
ax.set_xlabel("time [1/γ]") 
plt.tight_layout()
plt.savefig("two_continous_fockstate.svg") #hide
nothing #hide
```
![alt text](two_continous_fockstate.svg)





[^1]: [Heuck2020Photon-photonCavities](@cite)
[^2]: [Baragiola2012N-PhotonSystem](@cite)
