# WaveguideQED.jl
<a href="https://qojulia.github.io/WaveguideQED.jl/dev/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Documentation of latest stable version"></a> 
<a href="https://codecov.io/gh/qojulia/WaveguideQED.jl"><img src="https://img.shields.io/codecov/c/gh/qojulia/WaveguideQED.jl?label=codecov" alt="Test coverage from codecov"></a>


A Julia package for simulating quantum states of photon wavepackets using a discrete-time formalism [Phys. Rev. A 101, 042322](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.101.042322). The package works as an extension to [QuantumOptics.jl](https://qojulia.org/) where bases and operators from WaveguideQED.jl can be used together with operators and bases from QuantumOpics.jl. 

### Example of usage:
Define a waveguide basis, containing a two-photon wavepacket for a time interval from 0 to 20 with timesteps of 0.2:


```julia
using WaveguideQED
times = 0:0.1:10
bw =  WaveguideBasis(2,times)
```
Define waveguide creation and annihilation operators from this basis:

```julia
w = destroy(bw)
wd = create(bw)
```

Combine with QuantumOptics.jl operators:

```julia
using QuantumOptics
bc = FockBasis(2)
a = destroy(bc)
ad = create(bc)
wda = a ⊗ wd
adw = ad ⊗ w
```

Finally, we can define an initial two-photon Gaussian wavepacket state with view_twophoton and zero photons in the cavity, a Hamiltonian, and simulate the evolution:


```julia
ξfun(t1,t2,σ1,σ2,t0) = sqrt(2/σ1) * (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ1^2)*sqrt(2/σ2)*(log(2)/pi)^(1/4)*exp(-2*log(2)*(t2-t0)^2/σ2^2)
ψ_cw = twophoton(bw,ξfun,1,1,5)
psi = fockstate(bc,0) ⊗ ψ_cw
dt = times[2] - times[1]
H = im*sqrt(1/dt)*(adw-wda)
ψ = waveguide_evolution(times, psi, H)
```

Plotting the two-photon state is also simple:


```julia
ψ_double = TwoPhotonView(ψ);
using PyPlot
fig,ax = subplots(1,1,figsize=(9,4.5))
plot_twophoton!(ax,ψ_double,times)
```

![alt text](./Examples/two_photon_contour.jpg?raw=true)
