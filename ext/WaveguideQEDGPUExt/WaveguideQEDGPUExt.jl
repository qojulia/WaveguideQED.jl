module WaveguideQEDGPUExt

using WaveguideQED
import WaveguideQED: _is_destroy,_is_create,_is_identity,twophoton_index
using CUDA
using QuantumOpticsBase
using Strided
import QuantumOpticsBase: _tp_matmul_first!,_tp_matmul_last!
import LinearAlgebra: mul!


include("gpu_waveguideoperator.jl")
include("gpu_cavitywaveguideoperator.jl")
include("gpu_waveguideinteraction.jl")

end
