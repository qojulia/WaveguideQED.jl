# Unreleased
- `waveguide_evolution` now propagates each time bin with the matrix exponential (Krylov: Lanczos, or Arnoldi for non-Hermitian `H`) instead of an ODE solver (keywords `tol`, `maxiter`, `tstops`, `substeps`, `hermitian`). The Hamiltonian is constant within a bin, so this is exact up to `tol` and faster on all documentation examples (e.g. the GPU example: 5 s instead of 40 s).
- ODE-solver keywords such as `alg`, `abstol`, and `reltol` are ignored with a warning.
- `waveguide_montecarlo` moved to a package extension: it is available once `using QuantumOptics` has been run.
- Time-dependent Hamiltonians (e.g. `TimeDependentSum`) have their time set with `set_time!` during `waveguide_evolution`; `tstops` marks coefficient jumps inside bins.
- `fast_unitary` is deprecated in favor of `waveguide_evolution`, which applies the exponential of each bin instead of a truncated Taylor series.
- `waveguide_evolution` is precompiled: the first call in a new session takes 0.2 s instead of 4.6 s.
- `WaveguideBasis(Np, times)` no longer warns about the `lengths` keyword when it is not given.

# v0.2.14 - 13. May. 2026
- Improved compatibility with ForwardDiff-based workflows by making waveguide time-bin indexing robust to Dual-valued solver times. The discrete waveguide index now uses the primal time value internally, avoiding `round(Int, Dual, ...)` and `Float64(::Dual)` errors without adding a ForwardDiff dependency.

# v0.2.13 - 12. Apr. 2026
- Extended GPU operator capabilities to `WaveguideInteraction`.

# v0.2.12 - 10. Mar. 2026
- Added the possibility of having multiple waveguides of different lengths ([#68](https://github.com/qojulia/WaveguideQED.jl/pull/68)).

# v0.2.11 - 29. Nov. 2025
- Compat bump so that WaveguideQED works with the latest QuantumOptics ([#60](https://github.com/qojulia/WaveguideQED.jl/pull/60)).

# v0.2.10 - 23. May. 2025
- Small bugfix when constructing tensor products between waveguide operators and `SpinBasis` operators.

# v0.2.9 - 28. Apr. 2025
- Added `NLevelWaveguideOperator`s that speed up tensor products of waveguide operators and `NLevelBasis` (transition) operators ([#58](https://github.com/qojulia/WaveguideQED.jl/pull/58)).
- GPU support ([#55](https://github.com/qojulia/WaveguideQED.jl/pull/55)).

# v0.2.8 - 30. Jan. 2025
- Added `expect_waveguide()` to get expectation values of waveguide operators.
- Small bugfix for multiplying `LazyTensor`s with waveguide operators in them.

# v0.2.7 - 18. Aug. 2024
- Bugfix for `WaveguideInteraction` operators using delayed operators and looping ([#51](https://github.com/qojulia/WaveguideQED.jl/pull/51)).

# v0.2.6 - 31. Jul. 2024
- Waveguides now loop, so that if the simulation time extends the length of the waveguide state, the time index starts over. This also works with delayed operators, allowing for the simulation of separated emitters, etc. See also [#49](https://github.com/qojulia/WaveguideQED.jl/issues/49).
- Documentation added that showcases the usage of the above loop mechanic. See: https://qojulia.github.io/WaveguideQED.jl/dev/time_delay/#emitters

# v0.2.5 - 5. Jun. 2024
- Waveguide states are no longer, by default, initialized as normalized states. This is to ensure compatibility with creating superpositions / beamsplitter operations. In relation to this, the twophoton initializer was adjusted so that two photons in the same waveguide now require a prefactor of 1/sqrt(2) in order to be normalized. Documentation was also updated. See [#46](https://github.com/qojulia/WaveguideQED.jl/issues/46) for details on normalization.
- Precompilation updated to be faster.

# v0.2.4 - 28. May. 2024
- `expect_waveguide` to calculate expectation values of waveguide operators.
- Updated how twophoton operators are initialized to allow for non-symmetric wavefunctions.
- Added documentation showing how to simulate a twophoton state.

# v0.2.3 - 20. Feb. 2024
- `expect(O,rho)` now returns a sum over all timebins, and waveguide operators can now be added together with `LazySum`. See also [#44](https://github.com/qojulia/WaveguideQED.jl/issues/44).

# v0.2.2 - 15. Dec. 2023
- Added compatibility with `TimeDependentOperators` in QuantumOptics.jl (see [QuantumOpticsBase.jl#104](https://github.com/qojulia/QuantumOpticsBase.jl/pull/104)).

# v0.2.1 - 30. Jun. 2023
- added `delay` property to WaveguideOperators which can be set with the related keyword `delay` when created. They allow for non-markovian dynamics / feedback.
- added `effective_hamiltonian` and in general the file InputOutput.jl to deal with systems with more complex input-output relations
- Documentation overhaul
- Thesis source files added to repository
- Possible upstreamable extension of LazyOperators being multiplied added
