function mul!(result::Ket{B1,A}, a::WaveguideInteraction, b::Ket{B2,A}, alpha, beta) where {B1<:Basis,B2<:Basis,A<:CuArray}
    b_data = Base.ReshapedArray(b.data, QuantumOpticsBase._comp_size(basis(b)), ())
    result_data = Base.ReshapedArray(result.data, QuantumOpticsBase._comp_size(basis(result)), ())

    if a.loc == 1
        return QuantumOpticsBase._tp_matmul_first!(result_data, a, b_data, alpha * a.factor, beta)
    elseif a.loc == ndims(b_data)
        return QuantumOpticsBase._tp_matmul_last!(result_data, a, b_data, alpha * a.factor, beta)
    end
    QuantumOpticsBase._tp_matmul_mid!(result_data, a, loc, b_data, alpha * a.factor, beta)
      
    result
end


function QuantumOpticsBase._tp_matmul_mid!(result, a::WaveguideInteraction, loc::Integer, b, α::Number, β::Number)
    sz_b_1 = 1
    for i in 1:loc-1
        sz_b_1 *= size(b,i)
    end
    sz_b_3 = 1
    for i in loc+1:ndims(b)
        sz_b_3 *= size(b,i)
    end

    # TODO: Perhaps we should avoid reshaping here... should be possible to infer
    # contraction index tuple sizes
    br = Base.ReshapedArray(b, (sz_b_1, size(b, loc), sz_b_3), ())
    result_r = Base.ReshapedArray(result, (sz_b_1, size(a, 1), sz_b_3), ())

    # Try to "minimize" the transpose for efficiency.
    move_left = sz_b_1 < sz_b_3
    perm = move_left ? (2,1,3) : (1,3,2)

    br_p = QuantumOpticsBase._tp_matmul_get_tmp(eltype(br), ((size(br, i) for i in perm)...,), :_tp_matmul_mid_in, br)
    @strided permutedims!(br_p, br, perm)
    #permutedims!(br_p, br, perm)

    result_r_p = QuantumOpticsBase._tp_matmul_get_tmp(eltype(result_r), ((size(result_r, i) for i in perm)...,), :_tp_matmul_mid_out, result_r)
    @strided permutedims!(result_r_p, result_r, perm)
    #β == 0.0 || permutedims!(result_r_p, result_r, perm)

    if move_left
        QuantumOpticsBase._tp_matmul_first!(result_r_p, a, br_p, α, β)
    else
        QuantumOpticsBase._tp_matmul_last!(result_r_p, a, br_p, α, β)
    end

    @strided permutedims!(result_r, result_r_p, perm)
    #permutedims!(result_r, result_r_p, perm)

    result
end



function QuantumOpticsBase.:_tp_matmul_first!(result::Base.ReshapedArray{T,N1,AR1,<:Tuple}, a::WaveguideInteraction, b::Base.ReshapedArray{T,N2,AR2,<:Tuple}, α::Number, β::Number) where {
    T,
    N1,
    AR1<:CuArray,
    N2,
    AR2<:CuArray,
}
    d_first = size(b, 1)
    d_rest = length(b)÷d_first
    bp = b.parent
    rp = result.parent
    br = reshape(bp, (d_first, d_rest))
    result_r = reshape(rp, (size(a, 1), d_rest))
    apply_first_op_gpu!(result_r,a.op1,a.op2,br,α,β)
    result
end

#Same as _tp_matmul_first! But indexed in another way.
function QuantumOpticsBase.:_tp_matmul_last!(result::Base.ReshapedArray{T,N1,AR1,<:Tuple}, a::WaveguideInteraction, b::Base.ReshapedArray{T,N2,AR2,<:Tuple}, α::Number, β::Number) where {
    T,
    N1,
    AR1<:CuArray,
    N2,
    AR2<:CuArray,
}
    d_last = size(b, ndims(b))
    d_rest = length(b)÷d_last
    bp = b.parent
    rp = result.parent
    br = reshape(bp, (d_rest, d_last))
    result_r = reshape(rp, (d_rest, size(a, 1)))
    apply_last_op_gpu!(result_r,a.op1,a.op2,br,α,β)
    result
end


function apply_last_op_gpu!(
    result,
    a::WaveguideCreate{B1,B2,1,idx1},
    bop::WaveguideDestroy{B1,B2,1,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}

    # Figure out shape:
    N, M = size(result, 1), size(result, 2)

    # read waveguide info
    nsteps       = a.basis_l.nsteps
    timeindex_a  = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b  = (bop.timeindex + bop.delay -1) % nsteps + 1

    factor   = a.factor*bop.factor
    
    # Decide thread-block shape
    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_destroy_1photon_last!(
        result, b, α, β,
        N, M,
        nsteps,
        timeindex_a,
        timeindex_b,
        factor,
        idx1, idx2
    )
    return result
end


function gpu_interaction_create_destroy_1photon_last!(
    result, b, α, β,
    N, M,            # matrix shape
    nsteps,
    timeindex_a,
    timeindex_b,
    factor,
    idx1, idx2
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    # (A) multiply result by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # (B) do the "create x destroy" increment only for col == 1
    if col == 1 + (idx1-1)*nsteps + timeindex_a
        @inbounds result[row, col] += α * factor * b[row, 1 + (idx2-1)*nsteps + timeindex_b]
    end
end


#STILL IN PROGRESS

function apply_last_op_gpu!(
    result,
    a::WaveguideCreate{B1,B2,2,idx1},
    bop::WaveguideDestroy{B1,B2,2,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}

    # Figure out shape:
    N, M = size(result, 1), size(result, 2)

    # read waveguide info
    nsteps       = a.basis_l.nsteps
    timeindex_a  = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b  = (bop.timeindex + bop.delay -1) % nsteps + 1

    factor   = a.factor*bop.factor
    
    # Decide thread-block shape
    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_destroy_2photon_last!(
        result, b, α, β,
        N, M,
        nsteps,
        timeindex_a,
        timeindex_b,
        factor,
        idx1, idx2
    )
    return result
end

function gpu_interaction_create_destroy_2photon_last!(
    result, b, α, β,
    N, M,
    nsteps,
    timeindex_a,
    timeindex_b,
    factor,
    idx1, idx2
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    # (A) multiply by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # Direct implementation of the CPU's waveguide_interaction_mul! for WaveguideCreate+WaveguideDestroy
    
    # 1. Single-photon component (direct correlation between a photon in create and destroy)
    # This corresponds to: result[1+(idx1-1)*nsteps+timeindex_a] += alpha*factor*b[1+(idx2-1)*nsteps+timeindex_b]
    if col == 1 + (idx1-1)*nsteps + timeindex_a
        @inbounds result[row, col] += α * factor * b[row, 1 + (idx2-1)*nsteps + timeindex_b]
    end
    
    # 2. Two-photon operations
    
    # First case: Create from TwoPhotonTimestepView in B to TwoWaveguideTimestepView in Result
    # This corresponds to:
    # twophotonview_b = TwoPhotonTimestepView(b,timeindex_b,nsteps,Nw*nsteps+1+(idx2-1)*(nsteps*(nsteps+1))÷2)
    # i,j = min(idx1,idx2),max(idx1,idx2)
    # index = (i-1)*Nw + j - (i*(i+1))÷2
    # twophotonview_r = TwoWaveguideTimestepView(result,timeindex_a,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx1)
    
    twophoton_b_offset = Nw*nsteps + 1 + (idx2-1)*(nsteps*(nsteps+1))÷2
    i, j = min(idx1, idx2), max(idx1, idx2)
    index = (i-1)*Nw + j - (i*(i+1))÷2
    mixed_r_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
    
    # If we're in the mixed waveguide section of the output
    if col >= mixed_r_offset && col < mixed_r_offset + nsteps*nsteps
        # This is roughly equivalent to the TwoWaveguideTimestepView in the CPU implementation
        col_offset = col - mixed_r_offset
        t1 = (col_offset ÷ nsteps) + 1
        t2 = (col_offset % nsteps) + 1
        
        # We're in the right time slot for timeindex_a 
        # and mapping from a two-photon state in a single waveguide
        if (idx1 == i && t1 == timeindex_a) || (idx1 != i && t2 == timeindex_a)
            # Figure out which time corresponds to the other photon
            other_time = (idx1 == i) ? t2 : t1
            
            # Get the index in the source vector
            twophoton_idx = twophoton_b_offset + twophoton_index_gpu(min(other_time, timeindex_b), nsteps, max(other_time, timeindex_b))
            
            @inbounds result[row, col] += α * factor * b[row, twophoton_idx]
        end
    end
    
    # Second case: Create from TwoWaveguideTimestepView in B to TwoPhotonTimestepView in Result
    # This corresponds to:
    # twophotonview_b = TwoWaveguideTimestepView(b,timeindex_b,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx2)
    # twophotonview_r = TwoPhotonTimestepView(result,timeindex_a,nsteps,1+Nw*nsteps+(idx1-1)*(nsteps*(nsteps+1))÷2)
    
    mixed_b_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
    twophoton_r_offset = 1 + Nw*nsteps + (idx1-1)*(nsteps*(nsteps+1))÷2
    
    # If we're in the two-photon section of the output for waveguide idx1
    if col >= twophoton_r_offset && col < twophoton_r_offset + (nsteps*(nsteps+1))÷2
        # This is roughly equivalent to the TwoPhotonTimestepView in the CPU implementation
        col_offset = col - twophoton_r_offset
        
        # We need to determine the two times this represents
        is_int, t1 = reverse_idx_gpu(col_offset, nsteps, 0)
        
        if is_int
            # It's a valid time pair
            t2 = col_offset - twophoton_index_gpu(Int(t1), nsteps, 0) + Int(t1)
            
            # Mapping from a mixed waveguide state (photons in two different waveguides)
            mixed_idx = 0
            if idx2 == i # idx2 is the first waveguide in ordering
                mixed_idx = mixed_b_offset + (t1-1)*nsteps + t2
            else # idx2 is the second waveguide in ordering
                mixed_idx = mixed_b_offset + (t2-1)*nsteps + t1
            end
            
            # Assuming t1 is associated with timeindex_a (which might not be correct)
            if t1 == timeindex_a || t2 == timeindex_a
                @inbounds result[row, col] += α * factor * b[row, mixed_idx]
            end
        end
    end
    
    # Third case: Handle the loop for other waveguides
    # This corresponds to:
    # @simd for k in filter(x -> x != idx1 && x != idx2, 1:Nw)
    #   m,l = min(k,idx1),max(k,idx1)
    #   ... Rest of the loop handling other waveguide interactions
    
    # This is hard to vectorize efficiently in a GPU kernel, but we can approximate
    # by checking if the current column corresponds to any valid mixed waveguide state
    for k in 1:Nw
        if k == idx1 || k == idx2
            continue
        end
        
        # Calculate index for (k, idx1) pair
        m, l = min(k, idx1), max(k, idx1)
        index_r = (m-1)*Nw + l - (m*(m+1))÷2
        r_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index_r-1)*nsteps^2
        
        # Calculate index for (k, idx2) pair
        i, j = min(k, idx2), max(k, idx2)
        index_b = (i-1)*Nw + j - (i*(i+1))÷2
        b_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index_b-1)*nsteps^2
        
        # If we're in a mixed waveguide section for (k, idx1)
        if col >= r_offset && col < r_offset + nsteps*nsteps
            col_offset = col - r_offset
            t1 = (col_offset ÷ nsteps) + 1
            t2 = (col_offset % nsteps) + 1
            
            # We're in the right time slot for timeindex_a
            if (m == idx1 && t1 == timeindex_a) || (m != idx1 && t2 == timeindex_a)
                # Figure out which time corresponds to the other photon
                other_time = (m == idx1) ? t2 : t1
                
                # Get the source index in the b vector for the other waveguide pair
                b_time1 = 0
                b_time2 = 0
                
                if i == k # k is the first waveguide in the (k, idx2) ordering
                    b_time1 = other_time
                    b_time2 = timeindex_b
                else # k is the second waveguide in the (k, idx2) ordering
                    b_time1 = timeindex_b
                    b_time2 = other_time
                end
                
                mixed_b_idx = b_offset + (b_time1-1)*nsteps + b_time2
                
                @inbounds result[row, col] += α * factor * b[row, mixed_b_idx]
            end
        end
    end
    
    # Special case for vacuum state when operating on same waveguide/time
    if idx1 == idx2 && timeindex_a == timeindex_b && col == 1
        @inbounds result[row, col] += α * factor * b[row, 1]
    end
    
    return
end

function apply_first_op_gpu!(
    result,
    a::WaveguideCreate{B1,B2,2,idx1},
    bop::WaveguideDestroy{B1,B2,2,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}

    # Figure out shape:
    N, M = size(result, 1), size(result, 2)

    # read waveguide info
    nsteps       = a.basis_l.nsteps
    timeindex_a  = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b  = (bop.timeindex + bop.delay -1) % nsteps + 1
    Nw           = get_number_of_waveguides(a.basis_l)

    factor   = a.factor*bop.factor
    


function gpu_interaction_create_destroy_2photon_first!(
    result, b, α, β,
    N, M,
    nsteps,
    timeindex_a,
    timeindex_b,
    factor,
    idx1, idx2,
    Nw
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    # (A) multiply by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # Direct implementation of the CPU's waveguide_interaction_mul! for WaveguideCreate+WaveguideDestroy
    # But transposed for first-mode operations
    
    # 1. Single-photon component (direct correlation between a photon in create and destroy)
    if row == 1 + (idx1-1)*nsteps + timeindex_a
        @inbounds result[row, col] += α * factor * b[1 + (idx2-1)*nsteps + timeindex_b, col]
    end
    
    # 2. Two-photon operations 
    
    # First case: Create from TwoPhotonTimestepView in B to TwoWaveguideTimestepView in Result
    twophoton_b_offset = Nw*nsteps + 1 + (idx2-1)*(nsteps*(nsteps+1))÷2
    i, j = min(idx1, idx2), max(idx1, idx2)
    index = (i-1)*Nw + j - (i*(i+1))÷2
    mixed_r_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
    
    # If we're in the mixed waveguide section of the output (transposed)
    if row >= mixed_r_offset && row < mixed_r_offset + nsteps*nsteps
        # This is roughly equivalent to the TwoWaveguideTimestepView in the CPU implementation
        row_offset = row - mixed_r_offset
        t1 = (row_offset ÷ nsteps) + 1
        t2 = (row_offset % nsteps) + 1
        
        # We're in the right time slot for timeindex_a
        # and mapping from a two-photon state in a single waveguide
        if (idx1 == i && t1 == timeindex_a) || (idx1 != i && t2 == timeindex_a)
            # Figure out which time corresponds to the other photon
            other_time = (idx1 == i) ? t2 : t1
            
            # Get the index in the source vector
            twophoton_idx = twophoton_b_offset + twophoton_index_gpu(min(other_time, timeindex_b), nsteps, max(other_time, timeindex_b))
            
            @inbounds result[row, col] += α * factor * b[twophoton_idx, col]
        end
    end
    
    # Second case: Create from TwoWaveguideTimestepView in B to TwoPhotonTimestepView in Result
    mixed_b_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
    twophoton_r_offset = 1 + Nw*nsteps + (idx1-1)*(nsteps*(nsteps+1))÷2
    
    # If we're in the two-photon section of the output for waveguide idx1
    if row >= twophoton_r_offset && row < twophoton_r_offset + (nsteps*(nsteps+1))÷2
        # This is roughly equivalent to the TwoPhotonTimestepView in the CPU implementation
        row_offset = row - twophoton_r_offset
        
        # We need to determine the two times this represents
        is_int, t1 = reverse_idx_gpu(row_offset, nsteps, 0)
        
        if is_int
            # It's a valid time pair
            t2 = row_offset - twophoton_index_gpu(Int(t1), nsteps, 0) + Int(t1)
            
            # Mapping from a mixed waveguide state (photons in two different waveguides)
            mixed_idx = 0
            if idx2 == i # idx2 is the first waveguide in ordering
                mixed_idx = mixed_b_offset + (t1-1)*nsteps + t2
            else # idx2 is the second waveguide in ordering
                mixed_idx = mixed_b_offset + (t2-1)*nsteps + t1
            end
            
            # Assuming t1 is associated with timeindex_a (which might not be correct)
            if t1 == timeindex_a || t2 == timeindex_a
                @inbounds result[row, col] += α * factor * b[mixed_idx, col]
            end
        end
    end
    
    # Third case: Handle the loop for other waveguides
    for k in 1:Nw
        if k == idx1 || k == idx2
            continue
        end
        
        # Calculate index for (k, idx1) pair
        m, l = min(k, idx1), max(k, idx1)
        index_r = (m-1)*Nw + l - (m*(m+1))÷2
        r_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index_r-1)*nsteps^2
        
        # Calculate index for (k, idx2) pair
        i, j = min(k, idx2), max(k, idx2)
        index_b = (i-1)*Nw + j - (i*(i+1))÷2
        b_offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index_b-1)*nsteps^2
        
        # If we're in a mixed waveguide section for (k, idx1)
        if row >= r_offset && row < r_offset + nsteps*nsteps
            row_offset = row - r_offset
            t1 = (row_offset ÷ nsteps) + 1
            t2 = (row_offset % nsteps) + 1
            
            # We're in the right time slot for timeindex_a
            if (m == idx1 && t1 == timeindex_a) || (m != idx1 && t2 == timeindex_a)
                # Figure out which time corresponds to the other photon
                other_time = (m == idx1) ? t2 : t1
                
                # Get the source index in the b vector for the other waveguide pair
                b_time1 = 0
                b_time2 = 0
                
                if i == k # k is the first waveguide in the (k, idx2) ordering
                    b_time1 = other_time
                    b_time2 = timeindex_b
                else # k is the second waveguide in the (k, idx2) ordering
                    b_time1 = timeindex_b
                    b_time2 = other_time
                end
                
                mixed_b_idx = b_offset + (b_time1-1)*nsteps + b_time2
                
                @inbounds result[row, col] += α * factor * b[mixed_b_idx, col]
            end
        end
    end
    
    # Special case for vacuum state when operating on same waveguide/time
    if idx1 == idx2 && timeindex_a == timeindex_b && row == 1
        @inbounds result[row, col] += α * factor * b[1, col]
    end
    
    return
end


    # read waveguide info
    nsteps       = a.basis_l.nsteps
    timeindex_a  = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b  = (bop.timeindex + bop.delay -1) % nsteps + 1
    Nw           = get_number_of_waveguides(a.basis_l)

    factor   = a.factor*bop.factor
    
    # Decide thread-block shape
    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_destroy_destroy_2photon_last!(
        result, b, α, β,
        N, M,
        nsteps,
        timeindex_a,
        timeindex_b,
        factor,
        idx1, idx2,
        Nw
    )
    return result
end

function gpu_interaction_destroy_destroy_2photon_last!(
    result, b, α, β,
    N, M,
    nsteps,
    timeindex_a,
    timeindex_b,
    factor,
    idx1, idx2,
    Nw
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    # (A) multiply by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # Only contribute to vacuum state (col == 1)
    if col == 1
        if idx1 == idx2
            # Same waveguide
            t2 = max(timeindex_a, timeindex_b)
            t1 = min(timeindex_a, timeindex_b)
            factor_sqrt2 = (timeindex_a == timeindex_b) ? sqrt(2) : 1.0
            
            offset = Nw*nsteps + 1 + (idx1-1)*(nsteps*(nsteps+1)>>>1)
            two_photon_idx = twophoton_index_gpu(t1, nsteps, t2)
            
            @inbounds result[row, col] += factor_sqrt2 * α * factor * b[row, offset + two_photon_idx]
        else
            # Different waveguides
            i, j = min(idx1, idx2), max(idx1, idx2)
            index = (i-1)*Nw + j - (i*(i+1))÷2
            offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
            
            t1 = idx1 >= idx2 ? timeindex_a : timeindex_b
            t2 = idx1 >= idx2 ? timeindex_b : timeindex_a
            
            cross_idx = offset + (t1-1)*nsteps + t2
            
            @inbounds result[row, col] += α * factor * b[row, cross_idx]
        end
    end
    return
end

function apply_first_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,2,idx1},
    bop::WaveguideDestroy{B1,B2,2,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # Implementation similar to last op but transposed
    N, M = size(result, 1), size(result, 2)
    nsteps = a.basis_l.nsteps
    timeindex_a = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b = (bop.timeindex + bop.delay -1) % nsteps + 1
    Nw = get_number_of_waveguides(a.basis_l)
    factor = a.factor*bop.factor
    
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_destroy_destroy_2photon_first!(
        result, b, α, β, N, M, nsteps, timeindex_a, timeindex_b, factor, idx1, idx2, Nw
    )
    return result
end

function gpu_interaction_destroy_destroy_2photon_first!(
    result, b, α, β, N, M, nsteps, timeindex_a, timeindex_b, factor, idx1, idx2, Nw
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    if !isone(β)
        @inbounds result[row, col] *= β
    end

    if row == 1
        if idx1 == idx2
            t2 = max(timeindex_a, timeindex_b)
            t1 = min(timeindex_a, timeindex_b)
            factor_sqrt2 = (timeindex_a == timeindex_b) ? sqrt(2) : 1.0
            
            offset = Nw*nsteps + 1 + (idx1-1)*(nsteps*(nsteps+1)>>>1)
            two_photon_idx = twophoton_index_gpu(t1, nsteps, t2)
            
            @inbounds result[row, col] += factor_sqrt2 * α * factor * b[offset + two_photon_idx, col]
        else
            i, j = min(idx1, idx2), max(idx1, idx2)
            index = (i-1)*Nw + j - (i*(i+1))÷2
            offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
            
            t1 = idx1 >= idx2 ? timeindex_a : timeindex_b
            t2 = idx1 >= idx2 ? timeindex_b : timeindex_a
            
            cross_idx = offset + (t1-1)*nsteps + t2
            
            @inbounds result[row, col] += α * factor * b[cross_idx, col]
        end
    end
    return
end

# 2-photon implementations for Create-Create
function apply_last_op_gpu!(
    result,
    a::WaveguideCreate{B1,B2,2,idx1},
    bop::WaveguideCreate{B1,B2,2,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # Figure out shape:
    N, M = size(result, 1), size(result, 2)

    # read waveguide info
    nsteps       = a.basis_l.nsteps
    timeindex_a  = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b  = (bop.timeindex + bop.delay -1) % nsteps + 1
    Nw           = get_number_of_waveguides(a.basis_l)

    factor   = a.factor*bop.factor
    
    # Decide thread-block shape
    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_create_2photon_last!(
        result, b, α, β,
        N, M,
        nsteps,
        timeindex_a,
        timeindex_b,
        factor,
        idx1, idx2,
        Nw
    )
    return result
end

function gpu_interaction_create_create_2photon_last!(
    result, b, α, β,
    N, M,
    nsteps,
    timeindex_a,
    timeindex_b,
    factor,
    idx1, idx2,
    Nw
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    # (A) multiply by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # Contribution only from vacuum state (b[row, 1])
    if idx1 == idx2
        # Same waveguide
        t2 = max(timeindex_a, timeindex_b)
        t1 = min(timeindex_a, timeindex_b)
        factor_sqrt2 = (timeindex_a == timeindex_b) ? sqrt(2) : 1.0
        
        offset = Nw*nsteps + 1 + (idx1-1)*(nsteps*(nsteps+1)>>>1)
        two_photon_idx = twophoton_index_gpu(t1, nsteps, t2)
        
        if col == offset + two_photon_idx
            @inbounds result[row, col] += factor_sqrt2 * α * factor * b[row, 1]
        end
    else
        # Different waveguides
        i, j = min(idx1, idx2), max(idx1, idx2)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
        
        t1 = idx1 >= idx2 ? timeindex_a : timeindex_b
        t2 = idx1 >= idx2 ? timeindex_b : timeindex_a
        
        cross_idx = offset + (t1-1)*nsteps + t2
        
        if col == cross_idx
            @inbounds result[row, col] += α * factor * b[row, 1]
        end
    end
    return
end

function apply_first_op_gpu!(
    result,
    a::WaveguideCreate{B1,B2,2,idx1},
    bop::WaveguideCreate{B1,B2,2,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # Implementation similar to last op but transposed
    N, M = size(result, 1), size(result, 2)
    nsteps = a.basis_l.nsteps
    timeindex_a = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b = (bop.timeindex + bop.delay -1) % nsteps + 1
    Nw = get_number_of_waveguides(a.basis_l)
    factor = a.factor*bop.factor
    
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_create_2photon_first!(
        result, b, α, β, N, M, nsteps, timeindex_a, timeindex_b, factor, idx1, idx2, Nw
    )
    return result
end

function gpu_interaction_create_create_2photon_first!(
    result, b, α, β, N, M, nsteps, timeindex_a, timeindex_b, factor, idx1, idx2, Nw
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    if !isone(β)
        @inbounds result[row, col] *= β
    end

    if idx1 == idx2
        t2 = max(timeindex_a, timeindex_b)
        t1 = min(timeindex_a, timeindex_b)
        factor_sqrt2 = (timeindex_a == timeindex_b) ? sqrt(2) : 1.0
        
        offset = Nw*nsteps + 1 + (idx1-1)*(nsteps*(nsteps+1)>>>1)
        two_photon_idx = twophoton_index_gpu(t1, nsteps, t2)
        
        if row == offset + two_photon_idx
            @inbounds result[row, col] += factor_sqrt2 * α * factor * b[1, col]
        end
    else
        i, j = min(idx1, idx2), max(idx1, idx2)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
        
        t1 = idx1 >= idx2 ? timeindex_a : timeindex_b
        t2 = idx1 >= idx2 ? timeindex_b : timeindex_a
        
        cross_idx = offset + (t1-1)*nsteps + t2
        
        if row == cross_idx
            @inbounds result[row, col] += α * factor * b[1, col]
        end
    end
    return
end

# 2-photon implementations for Destroy-Create
function apply_last_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,2,idx1},
    bop::WaveguideCreate{B1,B2,2,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # Figure out shape:
    N, M = size(result, 1), size(result, 2)

    # read waveguide info
    nsteps       = a.basis_l.nsteps
    timeindex_a  = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b  = (bop.timeindex + bop.delay -1) % nsteps + 1
    Nw           = get_number_of_waveguides(a.basis_l)

    factor   = a.factor*bop.factor
    
    # Decide thread-block shape
    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_destroy_create_2photon_last!(
        result, b, α, β,
        N, M,
        nsteps,
        timeindex_a,
        timeindex_b,
        factor,
        idx1, idx2,
        Nw
    )
    return result
end

function gpu_interaction_destroy_create_2photon_last!(
    result, b, α, β,
    N, M,
    nsteps,
    timeindex_a,
    timeindex_b,
    factor,
    idx1, idx2,
    Nw
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    # (A) multiply by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # Handle same waveguide and same time bin case 
    if idx1 == idx2 && timeindex_a == timeindex_b
        # Vacuum state contribution
        if col == 1
            @inbounds result[row, col] += α * factor * b[row, 1]
        end
        
        # Add contributions to all waveguide states
        # This is the critical part that was missing
        
        # Check if we're in a waveguide state column (all columns after the first)
        for j in 1:Nw
            waveguide_offset = 1 + (j-1)*nsteps
            
            # Only operate on columns that are waveguide columns
            if waveguide_offset <= col <= waveguide_offset+nsteps-1
                t_idx = col - waveguide_offset + 1
                
                if j == idx1 && t_idx == timeindex_a
                    # Same waveguide at matching time bin, apply factor of 2
                    @inbounds result[row, col] += α * factor * 2 * b[row, col]
                else
                    # Other time bins or waveguides just get the normal contribution
                    @inbounds result[row, col] += α * factor * b[row, col]
                end
            end
        end
    else
        # Different waveguides or time bins
        if col == 1 + (idx2-1)*nsteps + timeindex_b
            @inbounds result[row, col] += α * factor * b[row, 1 + (idx1-1)*nsteps + timeindex_a]
        end
    end
    return
end

function apply_first_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,2,idx1},
    bop::WaveguideCreate{B1,B2,2,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    N, M = size(result, 1), size(result, 2)
    nsteps = a.basis_l.nsteps
    timeindex_a = (a.timeindex + a.delay -1) % nsteps + 1
    timeindex_b = (bop.timeindex + bop.delay -1) % nsteps + 1
    Nw = get_number_of_waveguides(a.basis_l)
    factor = a.factor*bop.factor
    
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_destroy_create_2photon_first!(
        result, b, α, β, N, M, nsteps, timeindex_a, timeindex_b, factor, idx1, idx2, Nw
    )
    return result
end

function gpu_interaction_destroy_create_2photon_first!(
    result, b, α, β, N, M, nsteps, timeindex_a, timeindex_b, factor, idx1, idx2, Nw
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # Handle same waveguide and same time bin case
    if idx1 == idx2 && timeindex_a == timeindex_b
        # Vacuum state contribution
        if row == 1
            @inbounds result[row, col] += α * factor * b[1, col]
        end
        
        # Add contributions to all waveguide states (transposed from last version)
        # This is the critical part that was missing
        
        # Check if we're in a waveguide state row (all rows after the first)
        for j in 1:Nw
            waveguide_offset = 1 + (j-1)*nsteps
            
            # Only operate on rows that are waveguide rows
            if waveguide_offset <= row <= waveguide_offset+nsteps-1
                t_idx = row - waveguide_offset + 1
                
                if j == idx1 && t_idx == timeindex_a
                    # Same waveguide at matching time bin, apply factor of 2
                    @inbounds result[row, col] += α * factor * 2 * b[row, col]
                else
                    # Other time bins or waveguides just get the normal contribution
                    @inbounds result[row, col] += α * factor * b[row, col]
                end
            end
        end
    else
        # Different waveguides or time bins
        if row == 1 + (idx2-1)*nsteps + timeindex_b
            @inbounds result[row, col] += α * factor * b[1 + (idx1-1)*nsteps + timeindex_a, col]
        end
    end
    return
end
