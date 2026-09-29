# Unreleased
- `waveguide_evolution` now propagates each time bin with the matrix exponential (Krylov: Lanczos, or Arnoldi for non-Hermitian `H`) instead of an ODE solver (keywords `tol`, `maxiter`, `tstops`, `substeps`, `hermitian`). The Hamiltonian is constant within a bin, so this is exact up to `tol` and faster on all documentation examples (e.g. the GPU example: 5 s instead of 40 s).
- **Breaking:** the ODE solver has been removed, and with it the dependency on QuantumOptics.jl (and thereby on OrdinaryDiffEq/SciMLBase); WaveguideQED now only builds on QuantumOpticsBase.jl. ODE-solver keywords such as `alg`, `abstol`, and `reltol` are ignored with a warning.
- `waveguide_montecarlo` moved to a package extension: it is available once `using QuantumOptics` has been run.
- Time-dependent Hamiltonians (e.g. `TimeDependentSum`) have their time set with `set_time!` during `waveguide_evolution`; `tstops` marks coefficient jumps inside bins.

# v0.2.1 - 30. Jun. 2023
- added `delay` property to WaveguideOperators which can be set with the related keyword `delay` when created. They allow for non-markovian dynamics / feedback.
- added `effective_hamiltonian` and in general the file InputOutput.jl to deal with systems with more complex input-output relations
- Documentation overhaul
- Thesis source files added to repository
- Possible upstreamable extension of LazyOperators being multiplied added
