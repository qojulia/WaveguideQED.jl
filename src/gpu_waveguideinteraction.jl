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
    N3,
    AR1<:CuArray{T, N3, CUDA.DeviceMemory},
    N2,
    N4,
    AR2<:CuArray{T, N4, CUDA.DeviceMemory},
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
    N3,
    AR1<:CuArray{T, N3, CUDA.DeviceMemory},
    N2,
    N4,
    AR2<:CuArray{T, N4, CUDA.DeviceMemory},
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




function gpu_interaction_create_destroy_2photon_last!(
    result, b, α, β,
    N, M,
    nsteps,
    timeindex_a, timeindex_b,
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

    offset   = Nw*nsteps + 1 + (idx1-1)*(nsteps*(nsteps+1)>>>1)    
    is_int,r_idx = reverse_idx(col-offset,nsteps,timeindex_a)

    offset_b   = Nw*nsteps + 1 + (idx2-1)*(nsteps*(nsteps+1)>>>1)    
    
    # (B) single-photon part: if col == 1
    if col ==  1 + (idx1-1)*nsteps + timeindex_a
        @inbounds result[row,col] += α * factor_a * factor_b * b[row, 1 + (idx2-1)*nsteps + timeindex_b]
    
    elseif col-offset == twophoton_index(timeindex_a, nsteps, timeindex_a)
        colprime = col - offset - twophoton_index(timeindex_b, nsteps, 0)
        colb = colprime + startcol - 1
        @inbounds result[row, col] += sqrt(2)*α*factor * b[row, colb]

    elseif twophoton_index(timeindex, nsteps, timeindex)+1 <= col-offset <= twophoton_index(timeindex, nsteps, nsteps)
        colprime = col -offset - twophoton_index(timeindex, nsteps, 0)
        colb = colprime + startcol - 1
        
        @inbounds result[row, col] += α*factor * b[row, colb]

    elseif is_int
        colb = r_idx + startcol - 1
        @inbounds result[row, col] += α*factor * b[row, colb]
    else
    end
    


    return
end


