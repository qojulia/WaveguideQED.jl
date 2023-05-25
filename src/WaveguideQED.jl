module WaveguideQED

using QuantumOptics
using Strided
using UnsafeArrays
import LinearAlgebra: axpy!, dot, mul!, rmul!
import QuantumOpticsBase: create, dagger, destroy, expect, identityoperator, tensor

export TwoPhotonTimestepView,TwoWaveguideTimestepView,OnePhotonView,TwoPhotonView,TwoWaveguideView,
    WaveguideBasis,zerophoton,onephoton,twophoton,view_waveguide,get_waveguidetimeindex,set_waveguidetimeindex!,get_dt,get_nsteps,get_waveguide_location,get_waveguide_basis,get_number_of_waveguides,get_waveguide_operators,
    WaveguideOperator,WaveguideDestroy,WaveguideCreate,
    CavityWaveguideAbsorption,CavityWaveguideEmission,emission,absorption,
    WaveguideInteraction,
    waveguide_evolution,waveguide_montecarlo,CavityWaveguideOperator,
    detect_single_click,detect_single_click!,LazyTensorKet,LazyTensorBra,LazySumKet,get_all_projectors,detect_double_click,detect_double_click!,Detector,
    plot_twophoton!,
    destroy,create,tensor,⊗,dagger,identityoperator,expect

include("view.jl")
include("basis.jl")
include("WaveguideOperator.jl")
include("CavityWaveguideOperator.jl")
include("WaveguideInteraction.jl")
include("solver.jl")
include("should_upstream.jl")
include("detection.jl")
include("plotting.jl")
include("precompile.jl")
end
