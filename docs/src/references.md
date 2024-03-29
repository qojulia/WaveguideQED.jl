# Suggested readings

### Theory and background
The theory of time-bin continuous fockstates is introduced in [Heuck2020Photon-photonCavities](@cite) where it is used to derive the equations of motion for Waveguide QED systems. The numerical method in this library is heavily based on this approach, where instead of deriving the equations, the systems are solved by applying the operators themselves. The time-bin method is also used in [Heuck2020](@cite) and [Krastanov2022](@cite).

### Similar approaches

The following is a list of approaches that are trying to solve problems that can also be treated with this library.
* [HughesWQEDMPP2021](@cite) considers feedback in waveguide systems and uses a space-discretized waveguide picture with Monte Carlo trajectories
* [Fischer2018](@cite) relates many approaches to solving WaveguideQED problems with each other and also introduces a framework that aims to deal with similar problems through master equations, photon counting and tracing out the waveguide.
* The SLH formalism introduced in [Kiilerich2019](@cite) and [Kiilerich2019](@cite) uses cascaded cavities to simulate quantum pulses. Further work also includes: [Yang2022](@cite) [Christiansen2023](@cite)

### Papers where we reproduce results from
* The theoretical results in [DynamicalPhotonLodahl2022](@cite) are reproduced in [Scattering on a two-level system](@ref lodahl).
* The semi-infinite waveguide with feedback considered in is considered in [HughesWQEDMPP2021](@cite) and we reproduce the results in [Waveguide and Mirror](@ref feedback)
* The fano transport studied in [Joanesarson2020](@cite) and input-output relations formulated in [Xu2016](@cite) are described and studied in [Input-Output Relations](@ref input)

# References
```@bibliography
```
