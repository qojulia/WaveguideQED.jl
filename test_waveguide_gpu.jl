using WaveguideQED
using CUDA
using QuantumOpticsBase
using QuantumOptics
using Test
using LinearAlgebra

# Helper functions to convert to GPU
function to_gpu(a::Operator)
    data = cu(map(ComplexF32, a.data))
    gpu_operator = Operator(a.basis_l, a.basis_r, data)
    return gpu_operator
end

function to_gpu(a::Ket)
    data = cu(map(ComplexF32, a.data))
    gpu_operator = Ket(a.basis, data)
    return gpu_operator
end

println("Setting up test for WaveguideInteraction with the same waveguide/timebin...")

# Setup for failing create-destroy and destroy-create tests
bc = FockBasis(3)  # Small cavity basis
times = 0:0.1:4    # Smaller time range
bw = WaveguideBasis(2, 2, times)  # Two waveguides, two photons

# Create test states
psi_one = fockstate(bc, 0) ⊗ onephoton(bw, 1, x->exp(-(x-2)^2))
psi_one_gpu = to_gpu(psi_one)

# Test Destroy-Create (same waveguide and timeindex)
println("Testing Destroy-Create in same waveguide & timeindex...")
w1c = create(bw, 1)
w1d = destroy(bw, 1)

# Make them operate at same timeindex 
w1c.timeindex = 20  # middle of the range
w1d.timeindex = 20  # same as create

# Create interaction operator
H = identityoperator(bc) ⊗ (w1d * w1c)

# Results (CPU)
psi_result = copy(psi_one)
QuantumOptics.mul!(psi_result, H, psi_one, 1.0, 0.0)

# Results (GPU)
psi_gpu_result = to_gpu(copy(psi_one))
QuantumOptics.mul!(psi_gpu_result, H, psi_one_gpu, 1.0, 0.0)

# Check output
norm_diff = norm(Array(psi_gpu_result.data) - psi_result.data)
println("Norm of difference (should be near 0): $norm_diff")

# Print values to compare
for i in 1:min(20, length(psi_result.data))
    if abs(psi_result.data[i]) > 1e-10 || abs(Array(psi_gpu_result.data)[i]) > 1e-10
        println("Index $i: CPU=$(psi_result.data[i]) GPU=$(Array(psi_gpu_result.data)[i])")
    end
end

# Test for approximate equality
if isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
    println("✅ Destroy-Create test PASSED")
else
    println("❌ Destroy-Create test FAILED")
end

# Test Create-Destroy (same waveguide and timeindex)
println("\nTesting Create-Destroy in same waveguide & timeindex...")

# Create interaction operator
H = identityoperator(bc) ⊗ (w1c * w1d)

# Results (CPU)
psi_result = copy(psi_one)
QuantumOptics.mul!(psi_result, H, psi_one, 1.0, 0.0)

# Results (GPU)
psi_gpu_result = to_gpu(copy(psi_one))
QuantumOptics.mul!(psi_gpu_result, H, psi_one_gpu, 1.0, 0.0)

# Check output
norm_diff = norm(Array(psi_gpu_result.data) - psi_result.data)
println("Norm of difference (should be near 0): $norm_diff")

# Print values to compare
for i in 1:min(20, length(psi_result.data))
    if abs(psi_result.data[i]) > 1e-10 || abs(Array(psi_gpu_result.data)[i]) > 1e-10
        println("Index $i: CPU=$(psi_result.data[i]) GPU=$(Array(psi_gpu_result.data)[i])")
    end
end

# Test for approximate equality
if isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
    println("✅ Create-Destroy test PASSED")
else
    println("❌ Create-Destroy test FAILED")
end

# Additional test for LazyTensor combinations
println("\nTesting LazyTensor combinations...")

bc = FockBasis(3)
times = 0:0.1:3
bw = WaveguideBasis(2, times)

# Create operators
a = destroy(bc)
ad = create(bc)
a_gpu = to_gpu(dense(a))
ad_gpu = to_gpu(dense(ad))

w1c = create(bw, 1)
w1d = destroy(bw, 1)

# Create test states
psi = fockstate(bc, 1) ⊗ onephoton(bw, 1, x->exp(-(x-1.5)^2))
psi_gpu = to_gpu(psi)

# Test LazyTensor with interaction
H_cpu = identityoperator(bc) ⊗ (w1c * w1d)

# Results (CPU)
psi_result = copy(psi)
QuantumOptics.mul!(psi_result, H_cpu, psi, 1.0, 0.0)

# Results (GPU) 
psi_gpu_result = to_gpu(copy(psi))
QuantumOptics.mul!(psi_gpu_result, H_cpu, psi_gpu, 1.0, 0.0)

# Check output
norm_diff = norm(Array(psi_gpu_result.data) - psi_result.data)
println("Norm of difference (should be near 0): $norm_diff")

# Test for approximate equality
if isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
    println("✅ LazyTensor test PASSED")
else
    println("❌ LazyTensor test FAILED")
end

println("\nAll tests complete")