# CavityWaveguide.jl

CavityWaveguide.jl is a package for simulating continous fockstates in waveguides. It expands on [`QuantumOptics.jl`](https://qojulia.org/) by adding custom basises, operators, and routines for doing detection. 

### Dev docs
Added functionalities:
* [`WaveguideBasis`](@ref) for representing the waveguide space and the related generator functions: [`zerophoton`](@ref), [`onephoton`](@ref), and [`twophoton`](@ref). Also see [`view_onephoton`](@ref) and [`view_twophoton`](@ref) for viewing the waveguide data for plotting.
* [`WaveguideOperator`](@ref) which are specialized operators allowing efficient annihilation and creation operators at each timebin in the waveguide. They are created by giving a basis to [`CavityWaveguide.destroy`](@ref) and [`CavityWaveguide.create`](@ref)
* Since the interaction between the waveguide timebin mode $k$ and cavity/emitter is always given as: $a^\dagger w_k - a w_k^\dagger$ there are specially optimized functions for doing these operations called [`CavityWaveguideOperator`](@ref) which are created using a fockbasis and a waveguide basis and the functions [`emission`](@ref) and [`absorption`](@ref).
* [`Detector`](@ref), [`LazyTensorKet`](@ref), and [`LazySumKet`](@ref) together with [`detect_single_click`](@ref) and [`detect_double_click`](@ref) allows one to do a beamsplitter interference and subsequent detection on photons comming from two waveguides. 



```@meta
DocTestSetup = quote
    using CavityWaveguide
end
```

