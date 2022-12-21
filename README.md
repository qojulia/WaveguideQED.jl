# CavityWaveguide.jl

A julia package for simulating quantum states of photon wavepackets using a discrete-time formalism [Phys. Rev. A 101, 042322](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.101.042322). Package works as an extension to [QuantumOptics.jl](https://qojulia.org/) where basis and operators from CavityWaveguide.jl can be used togehter with operators and basis from QuantumOpics.jl. See also /Examples for a single and twophoton pulse interacting with a cavity with a Kerr-non linearity. 


Define a waveguide basis, containing a two photon wavepacket for a time interval 0 to 20 with 0.2 timesteps:


```jldoctest
julia> bw =  WaveguideBasis(2,0:0.2:20)
WaveguideBasis{2}([10303], 10302, 0, 101, 1)
```

Define waveguide creation and annihilation operators from this basis:

```jldoctest
julia> w =  destroy(bw)
WaveguideDestroy(dim=10303x10303)
  basis: WaveguideBasis{2}([10303], 10302, 0, 101, 1)
julia> w =  create(bw)
WaveguideCreate(dim=10303x10303)
  basis: WaveguideBasis{2}([10303], 10302, 0, 101, 1)
```

Combine with QuantumOptics.jl operators:

```jldoctest
julia> using QuantumOptics.jl
julia> bc = FockBasis(2);
julia> a = destroy(bc);
julia> ad = create(bc);
julia> wda = a \otimes wd
LazyTensor(dim=7959x7959)
  basis: [Fock(cutoff=2) ⊗ WaveguideBasis{2}([2653], 2652, 0, 51, 51)]
  operators: 2
  indices: [1,2]
julia> adw = ad \otimes w
LazyTensor(dim=7959x7959)
  basis: [Fock(cutoff=2) ⊗ WaveguideBasis{2}([2653], 2652, 0, 51, 51)]
  operators: 2
  indices: [1,2]
```

Finally, we can define an initial twophoton gaussian wavepacket state with view_twophoton and zero photons in the cavity, an Hamiltonian, and simulate the evolution:


```jldoctest
julia> ξfun(t::Number,σ::Number,t0::Number) = 1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2)/sqrt(0.2820947917738782);
julia> ξvec=ξfun.(times,1,5)*transpose(ξfun.(times,1,5));
julia> ψ_cw = Ket(bw);
julia> tmp = view_twophoton(ψ_cw);
julia> tmp .= ξvec;
julia> psi = fockstate(bc,0) ⊗ ψ_cw;
julia> H = LazySum(im*sqrt(param.γ/dt)*adw,-im*sqrt(param.γ/dt)*wda);
julia> ψ = waveguide_evolution(param.times, psi, H);
```

Plotting the twophoton state is also simple:


```jldoctest
julia> ψ_double = view_twophoton(ψ);
julia> ψ_double = ψ_double + ψ_double' - Diagonal(ψ_double);

julia> using PyPlot
julia> pygui(true);
julia> fig,ax = subplots(1,1,figsize=(9,4.5))
julia> xgrid = repeat(param.times',length(param.times),1)
julia> ygrid = repeat(param.times,1,length(param.times))
julia> ax.contourf(xgrid,ygrid,ψ_double.*conj(ψ_double),100)
julia> ax.set_aspect("equal", "box")
julia> ax.set_ylabel("time [1/γ]")
julia> ax.set_xlabel("time [1/γ]")
julia> ax.set_title("δ = $(param.δ)")   
julia> plt.tight_layout()
```

Giving the following plot:
![alt text](./Examples/two_photon_contour.jpg?raw=true)
