# WaveguideQED.jl

WaveguideQED.jl is a package for simulating continuous fockstates in waveguides. It expands on [`QuantumOptics.jl`](https://qojulia.org/) by adding a custom basis, operators, and routines for doing detection. 

### Dev docs
Added functionalities:
* [`WaveguideBasis`](@ref) for representing the waveguide space and the related generator functions: [`zerophoton`](@ref), [`onephoton`](@ref), and [`twophoton`](@ref). Also see [`OnePhotonView`](@ref), [`TwoPhotonView`](@ref), and [`plot_twophoton!`](@ref) for viewing the waveguide data for plotting. Note that [`WaveguideBasis`](@ref) can contain multiple waveguides.
* [`WaveguideOperator`](@ref) which are specialized operators allowing efficient annihilation and creation operators at each time-bin in the waveguide. They are created by giving a basis to [`WaveguideQED.destroy`](@ref) and [`WaveguideQED.create`](@ref)
* Since the interaction between the waveguide time-bin mode $k$ and cavity/emitter is given as: $a^\dagger w_k - a w_k^\dagger$ there are specially optimized functions for doing these operations called [`CavityWaveguideOperator`](@ref) which are created using a fockbasis and a waveguide basis and the functions [`emission`](@ref) and [`absorption`](@ref).
* [`Detector`](@ref), [`LazyTensorKet`](@ref), and [`LazySumKet`](@ref), together with [`detect_single_click`](@ref) and [`detect_double_click`](@ref) allow one to do a beamsplitter interference and subsequent detection of photons coming from two waveguides. 

![alt text](./animations/firstgif.gif)

```@meta
DocTestSetup = quote
    using WaveguideQED
end
```

