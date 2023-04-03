# WaveguideQED.jl
<a href="https://mabuni1998.github.io/WaveguideQED.jl/dev/"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Documentation of latest stable version"></a> 
<a href="https://codecov.io/gh/mabuni1998/WaveguideQED.jl"><img src="https://img.shields.io/codecov/c/gh/mabuni1998/WaveguideQED.jl?label=codecov" alt="Test coverage from codecov"></a>


A julia package for simulating quantum states of photon wavepackets using a discrete-time formalism [Phys. Rev. A 101, 042322](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.101.042322). Package works as an extension to [QuantumOptics.jl](https://qojulia.org/) where basises and operators from WaveguideQED.jl can be used together with operators and basises from QuantumOpics.jl. 

### Example of usage:
Define a waveguide basis, containing a two photon wavepacket for a time interval 0 to 20 with 0.2 timesteps:


```jldoctest
julia> times = 0:0.1:20
julia> bw =  WaveguideBasis(2,times)
```
Define waveguide creation and annihilation operators from this basis:

```jldoctest
julia> w =  destroy(bw)
julia> w =  create(bw)
```

Combine with QuantumOptics.jl operators:

```jldoctest
julia> using QuantumOptics.jl
julia> bc = FockBasis(2);
julia> a = destroy(bc);
julia> ad = create(bc);
julia> wda = a \otimes wd
```

Finally, we can define an initial twophoton gaussian wavepacket state with view_twophoton and zero photons in the cavity, an Hamiltonian, and simulate the evolution:


```jldoctest
julia> ξfun(t1,t2,σ1,σ2,t0) = sqrt(2/σ1)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ1^2)*sqrt(2/σ2)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t2-t0)^2/σ2^2)
julia> ψ_cw = twophoton(bw,ξfun,times,1,1,5);
julia> psi = fockstate(bc,0) ⊗ ψ_cw;
julia> dt = times[2] - times[1];
julia> H = im*sqrt(1/dt)*(adw-wda);
julia> ψ = waveguide_evolution(param.times, psi, H);
```

Plotting the twophoton state is also simple:


```jldoctest
julia> ψ_double = view_twophoton(ψ);

julia> using PyPlot
julia> pygui(true);
julia> fig,ax = subplots(1,1,figsize=(9,4.5))
julia> xgrid = repeat(times',length(times),1)
julia> ygrid = repeat(times,1,length(times))
julia> ax.contourf(xgrid,ygrid,ψ_double.*conj(ψ_double),100)
julia> ax.set_aspect("equal", "box")
julia> ax.set_ylabel("time [1/γ]")
julia> ax.set_xlabel("time [1/γ]")
julia> ax.set_title("δ = 0")   
julia> plt.tight_layout()
```

Giving the following plot:
![alt text](./Examples/two_photon_contour.jpg?raw=true)
