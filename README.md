# CavityWaveguide.jl

A julia package for simulating quantum states of photon wavepackets using a discrete-time formalism [Phys. Rev. A 101, 042322](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.101.042322). Package works as an extension to [QuantumOptics.jl](https://qojulia.org/) where basis and operators from CavityWaveguide.jl can be used togehter with operators and basis from QuantumOpics.jl.


Define a waveguide basis, containing a two photon wavepacket for a time interval 0 to 20 with 0.2 timesteps:


```jldoctest
julia> bw =  WaveguideBasis(2,0:0.2:20)
WaveguideBasis{2}([10303], 10302, 0, 101, 1)
```


