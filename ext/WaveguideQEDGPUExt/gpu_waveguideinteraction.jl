function mul!(result::Ket{B1,A}, a::WaveguideInteraction{B1,B2}, b::Ket{B2,A}, alpha, beta) where {B1<:Basis,B2<:Basis,A<:CuArray}
    @views b_data = reshape(b.data, QuantumOpticsBase._comp_size(basis(b)))
    @views result_data = reshape(result.data, QuantumOpticsBase._comp_size(basis(result)))

    if a.loc == 1
        return QuantumOpticsBase._tp_matmul_first!(result_data, a, b_data, alpha * a.factor, beta)
    elseif a.loc == ndims(b_data)
        return QuantumOpticsBase._tp_matmul_last!(result_data, a, b_data, alpha * a.factor, beta)
    end
    QuantumOpticsBase._tp_matmul_mid!(result_data, a, a.loc, b_data, alpha * a.factor, beta)
      
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



function QuantumOpticsBase.:_tp_matmul_first!(
    result :: CuReshapedOrCuArray{T,N1},
    a      :: WaveguideInteraction,
    b      :: CuReshapedOrCuArray{T,N2},
    α      :: Number,
    β      :: Number
) where {T,N1,N2}
    d_first = size(b, 1)
    d_rest = length(b)÷d_first
    bp = parent_of(b)
    rp = parent_of(result)
    br = reshape(bp, (d_first, d_rest))
    result_r = reshape(rp, (size(a.op1, 1), d_rest))
    apply_first_op_gpu!(result_r,a.op1,a.op2,br,α,β)
    result
end

#Same as _tp_matmul_first! But indexed in another way.
function QuantumOpticsBase.:_tp_matmul_last!(
    result :: CuReshapedOrCuArray{T,N1},
    a      :: WaveguideInteraction,
    b      :: CuReshapedOrCuArray{T,N2},
    α      :: Number,
    β      :: Number
) where {T,N1,N2}
    d_last = size(b, ndims(b))
    d_rest = length(b)÷d_last
    bp = parent_of(b)
    rp = parent_of(result)
    br = reshape(bp, (d_rest, d_last))
    result_r = reshape(rp, (d_rest, size(a.op1, 1)))
    apply_last_op_gpu!(result_r,a.op1,a.op2,br,α,β)
    result
end

# Helper functions for GPU implementation
@inline function twophoton_index_gpu(i, nsteps, j)
    return (i-1)*nsteps + j - ((i-1)*i)÷2
end

@inline function reverse_idx_gpu(C::Int, nsteps::Int, timeindex::Int)
    out = nsteps+1/2-sqrt(4*nsteps^2 - 8*C-4*nsteps + 8*timeindex + 1)
    is_int = isinteger(out)
    return is_int, is_int ? Int(out) : out
end

@inline function reverse_idx_waveguide_gpu(C::Int, nsteps::Int, timeindex::Int)
    out = (C-timeindex)/nsteps+1
    is_int = isinteger(out) && (timeindex <= C <= (nsteps-1)*nsteps + timeindex)
    return is_int, is_int ? Int(out) : out
end

# Helper kernel for operations that only scale by β
function gpu_interaction_scale_only!(result, β, N, M)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    if !isone(β)
        @inbounds result[row, col] *= β
    end
    return
end

# One-photon operations
# Create-Destroy implementations (1-photon)
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

function apply_first_op_gpu!(
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
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_destroy_1photon_first!(
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
    return
end

function gpu_interaction_create_destroy_1photon_first!(
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

    # (B) do the "create x destroy" increment only for appropriate row
    if row == 1 + (idx1-1)*nsteps + timeindex_a
        @inbounds result[row, col] += α * factor * b[1 + (idx2-1)*nsteps + timeindex_b, col]
    end
    return
end

# Destroy-Destroy implementations (1-photon)
function apply_last_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,1,idx1},
    bop::WaveguideDestroy{B1,B2,1,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # For destroy-destroy (1-photon), we just scale by β as in CPU implementation
    N, M = size(result, 1), size(result, 2)
    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_scale_only!(
        result, β, N, M
    )
    return result
end

function apply_first_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,1,idx1},
    bop::WaveguideDestroy{B1,B2,1,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # For destroy-destroy (1-photon), we just scale by β as in CPU implementation
    N, M = size(result, 1), size(result, 2)
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_scale_only!(
        result, β, N, M
    )
    return result
end

# Create-Create implementations (1-photon)
function apply_last_op_gpu!(
    result,
    a::WaveguideCreate{B1,B2,1,idx1},
    bop::WaveguideCreate{B1,B2,1,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # For create-create (1-photon), we just scale by β as in CPU implementation
    N, M = size(result, 1), size(result, 2)
    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_scale_only!(
        result, β, N, M
    )
    return result
end

function apply_first_op_gpu!(
    result,
    a::WaveguideCreate{B1,B2,1,idx1},
    bop::WaveguideCreate{B1,B2,1,idx2},
    b,
    α,
    β
) where {B1,B2,idx1,idx2}
    # For create-create (1-photon), we just scale by β as in CPU implementation
    N, M = size(result, 1), size(result, 2)
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_scale_only!(
        result, β, N, M
    )
    return result
end

# Destroy-Create implementations (1-photon)
function apply_last_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,1,idx1},
    bop::WaveguideCreate{B1,B2,1,idx2},
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

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_destroy_create_1photon_last!(
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

function apply_first_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,1,idx1},
    bop::WaveguideCreate{B1,B2,1,idx2},
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
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_destroy_create_1photon_first!(
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

function gpu_interaction_destroy_create_1photon_last!(
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

    # Only add contribution if same waveguide and same time bin
    if idx1 == idx2 && timeindex_a == timeindex_b && col == 1
        @inbounds result[row, col] += α * factor * b[row, 1]
    end
    return
end

function gpu_interaction_destroy_create_1photon_first!(
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

    # Only add contribution if same waveguide and same time bin
    if idx1 == idx2 && timeindex_a == timeindex_b && row == 1
        @inbounds result[row, col] += α * factor * b[1, col]
    end
    return
end

# Two-photon operations
# Create-Destroy implementations (2-photon)
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
    Nw           = get_number_of_waveguides(a.basis_l)

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
        idx1, idx2,
        Nw
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

    # Single photon part
    if col == 1 + (idx1-1)*nsteps + timeindex_a
        @inbounds result[row, col] += α * factor * b[row, 1 + (idx2-1)*nsteps + timeindex_b]
    end

    # Two-photon same waveguide part
    if idx1 == idx2
        offset = Nw*nsteps + 1 + (idx1-1)*(nsteps*(nsteps+1)>>>1)
        
        # Handle photons in same waveguide
        twophoton_view_idx = twophoton_index_gpu(timeindex_a, nsteps, timeindex_b)
        if col >= offset && col < offset + (nsteps*(nsteps+1))÷2
            col_offset = col - offset
            
            # Adding support for same time bin (requires sqrt(2) factor)
            if timeindex_a == timeindex_b && col_offset == twophoton_index_gpu(timeindex_a, nsteps, timeindex_a)
                # This is the special case for same time bin
                @inbounds result[row, col] += sqrt(2) * α * factor * b[row, 1]
            else
                # Check if this index corresponds to one of our time bins
                is_valid = false
                
                # Two cases: our time bin is first or second in the pair
                if col_offset == twophoton_index_gpu(timeindex_a, nsteps, timeindex_b)
                    is_valid = true
                elseif col_offset == twophoton_index_gpu(timeindex_b, nsteps, timeindex_a)
                    is_valid = true
                end
                
                if is_valid
                    @inbounds result[row, col] += α * factor * b[row, 1]
                end
            end
        end
    else
        # Different waveguides
        i, j = min(idx1, idx2), max(idx1, idx2)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
        
        if idx1 == i  # First waveguide (row major ordering in array)
            target_idx = offset + (timeindex_a-1)*nsteps + timeindex_b
            if col == 1 + (idx1-1)*nsteps + timeindex_a
                @inbounds result[row, col] += α * factor * b[row, 1 + (idx2-1)*nsteps + timeindex_b]
            end
        else  # Second waveguide (column major ordering in array)
            target_idx = offset + (timeindex_b-1)*nsteps + timeindex_a
            if col == 1 + (idx1-1)*nsteps + timeindex_a
                @inbounds result[row, col] += α * factor * b[row, 1 + (idx2-1)*nsteps + timeindex_b]
            end
        end
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
    
    # Decide thread-block shape
    blockx, blocky = 32, 8
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_destroy_2photon_first!(
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

    # Single photon part
    if row == 1 + (idx1-1)*nsteps + timeindex_a
        @inbounds result[row, col] += α * factor * b[1 + (idx2-1)*nsteps + timeindex_b, col]
    end

    # Two-photon same waveguide part
    if idx1 == idx2
        offset = Nw*nsteps + 1 + (idx1-1)*(nsteps*(nsteps+1)>>>1)
        
        # Handle photons in same waveguide
        twophoton_view_idx = twophoton_index_gpu(timeindex_a, nsteps, timeindex_b)
        if row >= offset && row < offset + (nsteps*(nsteps+1))÷2
            row_offset = row - offset
            
            # Adding support for same time bin (requires sqrt(2) factor)
            if timeindex_a == timeindex_b && row_offset == twophoton_index_gpu(timeindex_a, nsteps, timeindex_a)
                # This is the special case for same time bin
                @inbounds result[row, col] += sqrt(2) * α * factor * b[1, col]
            else
                # Check if this index corresponds to one of our time bins
                is_valid = false
                
                # Two cases: our time bin is first or second in the pair
                if row_offset == twophoton_index_gpu(timeindex_a, nsteps, timeindex_b)
                    is_valid = true
                elseif row_offset == twophoton_index_gpu(timeindex_b, nsteps, timeindex_a)
                    is_valid = true
                end
                
                if is_valid
                    @inbounds result[row, col] += α * factor * b[1, col]
                end
            end
        end
        
    else
        # Different waveguides
        i, j = min(idx1, idx2), max(idx1, idx2)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        offset = 1 + Nw*nsteps + Nw*(nsteps*(nsteps+1))÷2 + (index-1)*nsteps^2
        
        if idx1 == i  # First waveguide (row major ordering in array)
            target_idx = offset + (timeindex_a-1)*nsteps + timeindex_b
            if row == 1 + (idx1-1)*nsteps + timeindex_a
                @inbounds result[row, col] += α * factor * b[1 + (idx2-1)*nsteps + timeindex_b, col]
            end
        else  # Second waveguide (column major ordering in array)
            target_idx = offset + (timeindex_b-1)*nsteps + timeindex_a
            if row == 1 + (idx1-1)*nsteps + timeindex_a
                @inbounds result[row, col] += α * factor * b[1 + (idx2-1)*nsteps + timeindex_b, col]
            end
        end
    end
    return
end

# 2-photon implementations for Destroy-Destroy
function apply_last_op_gpu!(
    result,
    a::WaveguideDestroy{B1,B2,2,idx1},
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

    # Only contribute to vacuum state at matching time bins
    if idx1 == idx2 && timeindex_a == timeindex_b
        if col == 1
            @inbounds result[row, col] += α * factor * b[row, 1]
            
            # Add for all waveguides (same as CPU implementation)
            if idx1 == idx2
                # For each waveguide, handle its contribution
                for j in 1:Nw
                    waveguide_offset = 1 + (j-1)*nsteps
                    
                    if j == idx1
                        # Same waveguide needs factor of 2 at exactly the time bin
                        if timeindex_a > 1 && timeindex_a < nsteps
                            # Contribute to all time bins except matching time bin
                            for t in 1:timeindex_a-1
                                @inbounds result[row, waveguide_offset+t-1] += b[row, waveguide_offset+t-1]
                            end
                            for t in timeindex_a+1:nsteps
                                @inbounds result[row, waveguide_offset+t-1] += b[row, waveguide_offset+t-1]
                            end
                            @inbounds result[row, waveguide_offset+timeindex_a-1] += 2 * b[row, waveguide_offset+timeindex_a-1]
                        end
                    else
                        # Other waveguides just get all time bins
                        for t in 1:nsteps
                            @inbounds result[row, waveguide_offset+t-1] += b[row, waveguide_offset+t-1]
                        end
                    end
                end
            end
        end
    else
        # Different time bins
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

    if idx1 == idx2 && timeindex_a == timeindex_b
        if row == 1
            @inbounds result[row, col] += α * factor * b[1, col]
            
            # For each waveguide, handle its contribution (transposed from last version)
            if idx1 == idx2
                for j in 1:Nw
                    waveguide_offset = 1 + (j-1)*nsteps
                    
                    if j == idx1
                        if timeindex_a > 1 && timeindex_a < nsteps
                            for t in 1:timeindex_a-1
                                @inbounds result[waveguide_offset+t-1, col] += b[waveguide_offset+t-1, col]
                            end
                            for t in timeindex_a+1:nsteps
                                @inbounds result[waveguide_offset+t-1, col] += b[waveguide_offset+t-1, col]
                            end
                            @inbounds result[waveguide_offset+timeindex_a-1, col] += 2 * b[waveguide_offset+timeindex_a-1, col]
                        end
                    else
                        for t in 1:nsteps
                            @inbounds result[waveguide_offset+t-1, col] += b[waveguide_offset+t-1, col]
                        end
                    end
                end
            end
        end
    else
        if row == 1 + (idx2-1)*nsteps + timeindex_b
            @inbounds result[row, col] += α * factor * b[1 + (idx1-1)*nsteps + timeindex_a, col]
        end
    end
    return
end