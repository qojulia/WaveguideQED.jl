# GPU Support
For particularly large system sizes there can be substantial computational gains by running the simulations on a Graphical Processing Unit (GPU). In WaveguideQED.jl, to run on a GPU, all operators and states (kets) should be "put" on the GPU, meaning they should be converted into CUDA arrays. The WaveguideOperators will then dispatch on custom written GPU kernels, which will significantly speed up the computational time for large systems. The GPU code is written using [https://github.com/JuliaGPU/CUDA.jl](CUDA.jl) and it should be installed before the following code will work, see [https://github.com/JuliaGPU/CUDA.jl](CUDA.jl) for details.  In the following, we provide an example of how to utilize a GPU.

To convert operators and kets to GPU arrays we define the following convinience functions:
```julia
using CUDA #Load CUDA

function to_gpu(a::Operator)
    data = cu(map(ComplexF32,a.data))
    gpu_operator = Operator(a.basis_l,a.basis_r,data)
    return gpu_operator
end

function to_gpu(a::Ket)
    data = cu(map(ComplexF32,a.data))
    gpu_ket = Ket(a.basis,data)
    return gpu_ket
end
```

We then define the system, its operators, and initial state:
```julia
dt = 0.01
times = 0:dt:10
bw = WaveguideBasis(2, times)
bc = FockBasis(2)


a = destroy(bc)
ad = create(bc)

w = destroy(bw)
wd = create(bw)

ξ₁(t1,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ^2)
ξ₂(t1,t2,σ1,σ2,t0) = ξ₁(t1,σ1,t0) * ξ₁(t2,σ2,t0) 
width = 1
t0 = 5
ψ = fockstate(bc,0) ⊗ twophoton(bw,1,ξ₂,width,width,t0)
```

WaveguideOperators does not need to be converted to a gpu, since they are just symbolic objects representing the action of the waveguide operator, but all system operators and all states need to be converted. In the following, we thus convert the operators and initial states to gpu operators:

```julia
a_gpu = to_gpu(a)
ad_gpu = to_gpu(ad)
ψ_gpu = to_gpu(ψ)
```

It is then simple to combine operators to run simulations on the gpu:

```julia
γ = 1
adw_gpu = ad_gpu ⊗ w
wda_gpu = a_gpu ⊗ wd
H_gpu =  sqrt(γ/dt)*(wda_gpu+ adw_gpu)

ψ_out = waveguide_evolution(times,ψ_gpu,H_gpu)
```

