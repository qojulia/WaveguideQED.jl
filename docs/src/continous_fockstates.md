# Continous Fock States

The single photon continuous fock state can be defined as:

$$\begin{equation*}
    \ket{\psi} = W^\dagger(\xi) \ket{0} = \int_{t_0}^{t_{end}} \mathrm{d}t \ \xi(t) w^\dagger(t) \ket{0}
\end{equation*}$$

here $W^\dagger(\xi)$ creates a photon with the wavefunction $\xi(t)$. $w^\dagger(t)$ is the creation operator for a photon at time $t$, and it obeys the commutation relation: $\comm{w(t)}{w(t')} = \delta(t-t')$. The probability of observing a photon at time $t$ is given by: $\bra{0} w(t) \ket{\psi} = |\xi(t)|^2$. The wavefunction $\xi(t)$ thus describes the temporal distribution of the photon.

The heart of the photon time binning is discretizing the continuous fock state into time bins of width $\Delta t$. It is then assumed that an emitter/cavity will only interact with the discretized waveguide state one timebin at a time. Corresponding to a spectrally flat interaction between the waveguide and emitter/cavity. We thus discretize the annihilation and creation operators by taking[^1]:

$$\begin{equation*}
    w(t_k) = w(k \Delta t) \rightarrow  \frac{w_k}{\sqrt{\Delta t}} \ \ \  \text{with} \ \left[ w_j, w_k^\dagger \right ] = \delta_{jk}
\end{equation*}$$

where $w_k$ is the descritized operator and the factor of $1/\sqrt{\Delta t}$ assures the commutator relation in the limit of $\Delta t \rightarrow 0$. This means that the single photon continuous fock state becomes:

$$\begin{equation*}
    \ket{\psi} = \frac{1}{\sqrt{2}} \int_{t_0}^{t_{end}} d t^{\prime} \int_{t_0}^{t_{end}} d t \ \xi(t) \xi\left(t^{\prime}\right) w^\dagger(t) w^\dagger\left(t^{\prime}\right)|0\rangle  \rightarrow 
\sum_{k=1}^N \sqrt{\Delta t} \xi(t_k) w_k^\dagger \ket{\emptyset}
\end{equation*}$$

In `CavityWaveguide.jl`, the timebins above are represented as elements in arrays corresponding to each timebin. Let`s say you want to represent a single photon contionous fock state that starts at $t=0$ and ends at $t=10$ with $\Delta t = 0.1$. This can be done be first creating waveguide basis capable defined on such a timeinterval:

```jldoctest
julia> times = 0:0.1:10
julia> bw = WaveguideBasis(1,times)
```

Notice that the input for WaveguideBasis is `1` and `times`. `1` denotes the maximum excitation number of fockstates (currently can only be 1 or 2) and `times` the time interval over which the continous fock state is defined. To define the continous fockstate we need to give a wavefuntion $\xi$. In the following we define a gaussian wavefunction located around $t=5$ with a width of $\sigma = 2$:

```jldoctest
julia> ξ(t,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2)
julia> σ,t0 = 2,5
julia> ψ = onephoton(bw,ξ,times,σ,t0)
```

This state can be visuallized by:

```jldoctest
julia> viewed_state = view_onephoton(ψ)
julia> fig,ax = subplots(1,1,figsize=(9,4.5))
julia> ax.plot(times,viewed_state,"r-")
julia> ax.set_xlabel("Time [a.u]")
julia> ax.set_ylabel(L"$\xi(t)$")
julia> plt.tight_layout()
```
![alt text](one_continous_fockstate.jpg)


We can straightforwardly extend the above definition to a two-photon continuous state as[^2]:

$$\begin{align*}
\frac{1}{\sqrt{2}}\left[W^\dagger(\xi)\right]^2|0\rangle = \frac{1}{\sqrt{2}} \int_{t_0}^{t_{end}} d t^{\prime} \int_{t_0}^{t_{end}} d t \ \xi(t) \xi\left(t^{\prime}\right) w^\dagger(t) w^\dagger\left(t^{\prime}\right)|0\rangle  
\end{align*}$$

The state is now defined over two times, which now describes the probability of observing photon A at time $t$ and photon B at time $t'$. In this case, the state is a product state, and both probabilities are described by the single photon wavefunction $\xi(t)$. 

For the timebinning we now have:

$$\begin{align*}
\frac{1}{\sqrt{2}}\left[W^\dagger(\xi)\right]^2|0\rangle &= \frac{1}{\sqrt{2}} \int_{t_0}^{t_{end}} d t^{\prime} \int_{t_0}^{t_{end}} d t \ \xi(t) \xi\left(t^{\prime}\right) w^\dagger(t) w^\dagger\left(t^{\prime}\right)|0\rangle \\
& \rightarrow \frac{1}{\sqrt{2}} \sum_{i=1}^N \sum_{k=1}^N \xi\left(t_i\right) \xi\left(t_k\right) w^\dagger\left(t_i\right) w^{\dagger}\left(t_k\right)|0\rangle \\
& =\frac{1}{\sqrt{2}} \sum_{i=1}^N \sum_{k \neq i}^N \xi\left(t_i\right) \xi\left(t_k\right) w^{\dagger}\left(t_i\right) w^{\dagger}\left(t_k\right)|0\rangle+\sum_{i=1}^N \xi\left(t_i\right) \xi\left(t_i\right)\left|2 t_i\right\rangle \\
& =\frac{2}{\sqrt{2}} \sum_{i=1}^N \sum_{k>i}^N \xi\left(t_i\right) \xi\left(t_k\right)\left|1_{t_i} 1_{t_k}\right\rangle+\sum_{i=1}^N \xi\left(t_i\right) \xi\left(t_i\right)\left|2 t_i\right\rangle \\
& =\sqrt{2} \sum_{i=1}^N \sum_{k > i}^N \xi\left(t_i\right) \xi\left(t_k\right) \mid 1_{t_i} 1_{t_k}\rangle + \sum_{i=1}^N \xi\left(t_i\right) \xi\left(t_i\right)\left|2 t_i\right\rangle
\end{align*}$$


where the sum is allowed to run over only half of the times due to the symmetry of the photons (it's equivalent having one photon at time bin k and then one photon at time bin j or one photon at time bin j and then one photon at time bin k). This is how the twophoton state is saved in the underlying arrays and we can define a twophoton basis as and corresponding twophoton state as:

```jldoctest
julia> bw = WaveguideBasis(2,times)
julia> ξ2(t1,t2,σ,t0) = ξ(t1,σ,t0)*ξ(t2,σ,t0)
julia> σ,t0 = 2,5
julia> ψ = twophoton(bw,ξ2,times,σ,t0)
```

Notice that we here defined a symmetric twophoton product state 

```
julia> viewed_state = view_onephoton(ψ)
julia> fig,ax = subplots(1,1,figsize=(9,4.5))
julia> ax.plot(times,viewed_state,"r-")
julia> ax.set_xlabel("Time [a.u]")
julia> ax.set_ylabel(L"$\xi(t)$")
julia> plt.tight_layout()
```
![alt text](two_continous_fockstate.jpg)


[^1]: [Heuck2020Photon-photonCavities](@cite)
[^2]: [Baragiola2012N-PhotonSystem](@cite)
