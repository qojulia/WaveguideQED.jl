module WaveguideQED

using QuantumOptics
using Strided
using UnsafeArrays
using FFTW
using Test
import LinearAlgebra: axpy!, dot, mul!, rmul!,I
import QuantumOpticsBase: create, dagger, destroy, expect, identityoperator, tensor,set_time!

export TwoPhotonTimestepView,TwoWaveguideTimestepView,OnePhotonView,TwoPhotonView,TwoWaveguideView,
    WaveguideBasis,zerophoton,onephoton,twophoton,view_waveguide,get_waveguidetimeindex,set_waveguidetimeindex!,get_dt,get_nsteps,get_waveguide_location,get_waveguide_basis,get_number_of_waveguides,get_waveguide_operators,expect_waveguide,
    WaveguideOperator,WaveguideDestroy,WaveguideCreate,
    CavityWaveguideAbsorption,CavityWaveguideEmission,emission,absorption,
    WaveguideInteraction,
    waveguide_evolution,waveguide_montecarlo,fast_unitary,
    CavityWaveguideOperator,
    detect_single_click,detect_single_click!,LazyTensorKet,LazyTensorBra,LazySumKet,get_all_projectors,detect_double_click,detect_double_click!,Detector,
    plot_twophoton!,
    WaveguideTransform,WaveguideTransformEvolution,effective_hamiltonian,fftket,
    destroy,create,tensor,⊗,dagger,identityoperator,set_time!,
    linear_index

include("view.jl")
include("basis.jl")
include("WaveguideOperator.jl")
include("CavityWaveguideOperator.jl")
include("NLevelWaveguideOperator.jl")
include("WaveguideInteraction.jl")
include("solver.jl")
include("detection.jl")
include("plotting.jl")
include("InputOutput.jl")
include("precompile.jl")
end
