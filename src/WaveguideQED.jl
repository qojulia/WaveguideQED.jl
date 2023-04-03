module WaveguideQED

using QuantumOptics,UnsafeArrays,Strided
import LinearAlgebra: mul!, rmul!,axpy!
import QuantumOptics.QuantumOpticsBase: dagger,tensor,create,destroy
import QuantumOptics: identityoperator

export TwoPhotonTimestepView,TwophotonView,TwoWaveguideTimestepView,OnePhotonView,TwoPhotonView,
    WaveguideBasis,zerophoton,onephoton,twophoton,view_waveguide,get_waveguidetimeindex,set_waveguidetimeindex!,get_nsteps,get_waveguide_location,get_waveguide_basis,get_waveguide_operators,
    WaveguideOperator,WaveguideDestroy,WaveguideCreate,
    CavityWaveguideAbsorption,CavityWaveguideEmission,emission,absorption,
    WaveguideInteraction,
    waveguide_evolution,waveguide_montecarlo,CavityWaveguideOperator,
    detect_single_click,detect_single_click!,LazyTensorKet,LazyTensorBra,LazySumKet,get_all_projectors,detect_double_click,detect_double_click!,Detector,
    plot_twophoton!,
    destroy,create,tensor,⊗,dagger,identityoperator

include("basis.jl")
include("view.jl")
include("WaveguideOperator.jl")
include("CavityWaveguideOperator.jl")
include("WaveguideInteraction.jl")
include("solver.jl")
include("should_upstream.jl")
include("detection.jl")
include("plotting.jl")
end
