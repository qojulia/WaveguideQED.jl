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
# Setup bases
bc = FockBasis(4)  # Small cavity basis for tests
times = 0:0.1:5    # Smaller time range for tests
bw = WaveguideBasis(1,2, times)  # Two waveguides

# Create operators
psi = fockstate(bc, 0) ⊗ onephoton(bw, 1, x->exp(-(x-2.5)^2))
psi_gpu = to_gpu(psi)
# Create waveguide operators
w1c = create(bw, 1)
w1d = destroy(bw, 1)
w2c = create(bw, 2)
w2d = destroy(bw, 2)

# Create interaction operator (CPU)
H = identityoperator(bc) ⊗ (w1c * w1d)

# Results (CPU)
psi_result = copy(psi)
QuantumOptics.mul!(psi_result, H, psi, 1.0, 0.0)

# Results (GPU)
psi_gpu_result = to_gpu(copy(psi));
QuantumOptics.mul!(psi_gpu_result, H, psi_gpu, 1.0, 0.0)

# Compare
@test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5

# Test cross-waveguide interaction
H_cross = identityoperator(bc) ⊗ ( w1c * w2d)

# Results (CPU)
psi_result = copy(psi)
QuantumOptics.mul!(psi_result, H_cross, psi, 1.0, 0.0)

# Results (GPU)
psi_gpu_result = to_gpu(copy(psi))
QuantumOptics.mul!(psi_gpu_result, H_cross, psi_gpu, 1.0, 0.0)

# Compare
@test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5

@testset "WaveguideQEDGPUExt - WaveguideInteraction" begin
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

    @testset "One-photon operators" begin
        # Setup bases
        bc = FockBasis(4)  # Small cavity basis for tests
        times = 0:0.1:5    # Smaller time range for tests
        bw = WaveguideBasis(1,2, times)  # Two waveguides
        
        # Create operators
        psi = fockstate(bc, 0) ⊗ onephoton(bw, 1, x->exp(-(x-2.5)^2))
        psi.data .= 1
        psi_gpu = to_gpu(psi)
        
        # Test different combinations
        @testset "Create-Destroy" begin
            # Create waveguide operators
            w1c = create(bw, 1)
            w1d = destroy(bw, 1)
            w2c = create(bw, 2)
            w2d = destroy(bw, 2)
            
            # Create interaction operator (CPU)
            H = identityoperator(bc) ⊗ ( w1c * w1d)
            # Results (CPU)
            psi_result = copy(psi)
            QuantumOptics.mul!(psi_result, H, psi, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi))
            QuantumOptics.mul!(psi_gpu_result, H, psi_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
            
            # Test cross-waveguide interaction
            H_cross = identityoperator(bc) ⊗ ( w1c * w2d)
            
            # Results (CPU)
            psi_result = copy(psi)
            QuantumOptics.mul!(psi_result, H_cross, psi, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
        end
        
        @testset "Destroy-Destroy" begin
            # Create waveguide operators
            w1d = destroy(bw, 1)
            w2d = destroy(bw, 2)
            
            # Create interaction operator
            H = identityoperator(bc) ⊗ ( w1d * w1d)
            
            # Results (CPU)
            psi_result = copy(psi)
            QuantumOptics.mul!(psi_result, H, psi, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi))
            QuantumOptics.mul!(psi_gpu_result, H, psi_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
            
            # Cross waveguide
            H_cross = identityoperator(bc) ⊗ ( w1d * w2d)
            
            # Results (CPU)
            psi_result = copy(psi)
            QuantumOptics.mul!(psi_result, H_cross, psi, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
        end
        
        @testset "Create-Create" begin
            # Create waveguide operators
            w1c = create(bw, 1)
            w2c = create(bw, 2)
            
            # Create state with vacuum component
            psi_vac = (fockstate(bc, 0) ⊗ zerophoton(bw))
            psi_vac_gpu = to_gpu(psi_vac)
            
            # Create interaction operator
            H = identityoperator(bc) ⊗ ( w1c * w1c)
            
            # Results (CPU)
            psi_result = copy(psi_vac)
            QuantumOptics.mul!(psi_result, H, psi_vac, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_vac))
            QuantumOptics.mul!(psi_gpu_result, H, psi_vac_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
            
            # Cross waveguide
            H_cross = identityoperator(bc) ⊗ ( w1c * w2c)
            
            # Results (CPU)
            psi_result = copy(psi_vac)
            QuantumOptics.mul!(psi_result, H_cross, psi_vac, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_vac))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_vac_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
        end
        
        @testset "Destroy-Create" begin
            # Create waveguide operators
            w1c = create(bw, 1)
            w1d = destroy(bw, 1)
            w2c = create(bw, 2)
            w2d = destroy(bw, 2)
            
            # Create interaction operator
            H = identityoperator(bc) ⊗ ( w1d * w1c)
            
            # Results (CPU)
            psi_result = copy(psi)
            QuantumOptics.mul!(psi_result, H, psi, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi))
            QuantumOptics.mul!(psi_gpu_result, H, psi_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
            
            # Cross waveguide
            H_cross = identityoperator(bc) ⊗ ( w1d * w2c)
            
            # Results (CPU)
            psi_result = copy(psi)
            QuantumOptics.mul!(psi_result, H_cross, psi, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_gpu, 1.0, 0.0)
            
            # Compare
            @test Array(psi_gpu_result.data) ≈ psi_result.data rtol=1e-5
        end
    end
    
    @testset "Two-photon operators" begin
        # Setup bases
        bc = FockBasis(3)  # Small cavity basis
        times = 0:0.1:4    # Smaller time range
        bw = WaveguideBasis(2, 2, times)  # Two waveguides, two photons
        
        # Create test states
        psi_vac = fockstate(bc, 0) ⊗ zerophoton(bw)
        psi_vac_gpu = to_gpu(psi_vac)
        
        psi_one = fockstate(bc, 0) ⊗ twophoton(bw, 1, (x,y)->exp(-(x-2)^2-(y-3)^2))
        psi_one = fockstate(bc, 0) ⊗ twophoton(bw, 1, (x,y)->exp(-(x-2)^2-(y-3)^2))
        
        #psi_one.data .= 1
        psi_one_gpu = to_gpu(psi_one)
        
        function gaussian_twophoton(t1, t2)
            return exp(-(t1-1.5)^2) * exp(-(t2-2.5)^2)
        end
        
        psi_two = fockstate(bc, 0) ⊗ twophoton(bw, 1, gaussian_twophoton)
        psi_two.data .= 1
        psi_two_gpu = to_gpu(psi_two)
        
        @testset "Create-Destroy" begin
            # Create waveguide operators
            w1c = create(bw, 1)  # 2-photon operators
            w1d = destroy(bw, 1)
            w2c = create(bw, 2)
            w2d = destroy(bw, 2)
            
            # Create interaction operator
            H = identityoperator(bc) ⊗ ( w1c * w1d)
            
            # Results (CPU)
            psi_result = copy(psi_two)
            QuantumOptics.mul!(psi_result, H, psi_two, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_two))
            QuantumOptics.mul!(psi_gpu_result, H, psi_two_gpu, 1.0, 0.0)
            
            println(findall(x->x!=0,abs.(Array(psi_gpu_result.data))))
            println(findall(x->x!=0,abs.(psi_result.data)))

            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
            
            # Cross waveguide
            H_cross = identityoperator(bc) ⊗ ( w1c * w2d)
            
            # Results (CPU)
            psi_result = copy(psi_two)
            QuantumOptics.mul!(psi_result, H_cross, psi_two, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_two))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_two_gpu, 1.0, 0.0)
            
            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
        end
        
        @testset "Destroy-Destroy" begin
            # Create waveguide operators
            w1d = destroy(bw, 1)
            w2d = destroy(bw, 2)
            
            # Create interaction operator
            H = identityoperator(bc) ⊗ ( w1d * w1d)
            
            # Results (CPU)
            psi_result = copy(psi_two)
            QuantumOptics.mul!(psi_result, H, psi_two, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_two))
            QuantumOptics.mul!(psi_gpu_result, H, psi_two_gpu, 1.0, 0.0)
            
            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
            
            # Cross waveguide
            H_cross = identityoperator(bc) ⊗ ( w1d * w2d)
            
            # Results (CPU)
            psi_result = copy(psi_two)
            QuantumOptics.mul!(psi_result, H_cross, psi_two, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_two))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_two_gpu, 1.0, 0.0)
            
            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
        end
        
        @testset "Create-Create" begin
            # Create waveguide operators
            w1c = create(bw, 1)
            w2c = create(bw, 2)
            
            # Create interaction operator
            H = identityoperator(bc) ⊗ ( w1c * w1c)
            
            # Results (CPU)
            psi_result = copy(psi_vac)
            QuantumOptics.mul!(psi_result, H, psi_vac, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_vac))
            QuantumOptics.mul!(psi_gpu_result, H, psi_vac_gpu, 1.0, 0.0)
            
            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
            
            # Cross waveguide
            H_cross = identityoperator(bc) ⊗ ( w1c * w2c)
            
            # Results (CPU)
            psi_result = copy(psi_vac)
            QuantumOptics.mul!(psi_result, H_cross, psi_vac, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_vac))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_vac_gpu, 1.0, 0.0)
            
            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
        end
        
        @testset "Destroy-Create" begin
            # Create waveguide operators
            w1c = create(bw, 1)
            w1d = destroy(bw, 1)
            w2c = create(bw, 2)
            w2d = destroy(bw, 2)
            
            # Create interaction operator
            H = identityoperator(bc) ⊗ ( w1d * w1c)
            
            # Results (CPU)
            psi_result = copy(psi_one)
            QuantumOptics.mul!(psi_result, H, psi_one, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_one))
            QuantumOptics.mul!(psi_gpu_result, H, psi_one_gpu, 1.0, 0.0)
            
            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
            
            # Cross waveguide
            H_cross = identityoperator(bc) ⊗ ( w1d * w2c)
            
            # Results (CPU)
            psi_result = copy(psi_one)
            QuantumOptics.mul!(psi_result, H_cross, psi_one, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi_one))
            QuantumOptics.mul!(psi_gpu_result, H_cross, psi_one_gpu, 1.0, 0.0)
            
            # Compare
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
        end
    end

    @testset "Complex operators and LazyTensor combinations" begin
        # Setup bases
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
        H_gpu = identityoperator(bc) ⊗ (w1c * w1d)
        QuantumOptics.mul!(psi_gpu_result, H_gpu, psi_gpu, 1.0, 0.0)
        
        # Compare
        @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
        
        # Test sum of interactions
        H_sum_cpu = identityoperator(bc) ⊗ ((w1c * w1d) + 0.5 * (w1d * w1c))
        
        # Results (CPU)
        psi_result = copy(psi)
        QuantumOptics.mul!(psi_result, H_sum_cpu, psi, 1.0, 0.0)
        
        # Results (GPU)
        psi_gpu_result = to_gpu(copy(psi))
        QuantumOptics.mul!(psi_gpu_result, H_sum_cpu, psi_gpu, 1.0, 0.0)
        
        # Compare
        @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
    end
    
    @testset "Performance metrics" begin
        # Only run timing tests if not in CI environment
        if !haskey(ENV, "CI")
            # Setup larger bases for meaningful performance comparison
            bc = FockBasis(8)
            times = 0:0.05:10
            bw = WaveguideBasis(2, times)
            
            # Create operators
            w1c = create(bw, 1)
            w1d = destroy(bw, 1)
            H = identityoperator(bc) ⊗ ( w1c * w1d)
            
            # Create test state
            psi = fockstate(bc, 1) ⊗ onephoton(bw, 1, x->exp(-(x-5)^2))
            psi.data .= 1
            psi_gpu = to_gpu(psi)
            
            # Results (CPU)
            psi_result = copy(psi)
            t_cpu = @elapsed QuantumOptics.mul!(psi_result, H, psi, 1.0, 0.0)
            
            # Results (GPU)
            psi_gpu_result = to_gpu(copy(psi))
            # Warmup
            QuantumOptics.mul!(psi_gpu_result, H, psi_gpu, 1.0, 0.0)
            CUDA.synchronize()
            
            t_gpu = @elapsed begin
                QuantumOptics.mul!(psi_gpu_result, H, psi_gpu, 1.0, 0.0)
                CUDA.synchronize()
            end
            
            println("CPU time: $t_cpu seconds")
            println("GPU time: $t_gpu seconds")
            println("Speedup: $(t_cpu/t_gpu)x")
            
            # We don't assert on speedup as it depends on hardware
            # but we verify correctness
            @test isapprox(Array(psi_gpu_result.data), psi_result.data; rtol=1e-5, atol=1e-5)
        end
    end
end