const InteractionWaveguideOp = Union{WaveguideCreate,WaveguideDestroy}


function mul!(result::Ket{B1,A}, a::WaveguideInteraction{B1,B2}, b::Ket{B2,A}, alpha, beta) where {B1<:Basis,B2<:Basis,A<:CuArray}
    b_data = Base.ReshapedArray(b.data, QuantumOpticsBase._comp_size(basis(b)), ())
    result_data = Base.ReshapedArray(result.data, QuantumOpticsBase._comp_size(basis(result)), ())

    if a.loc == 1
        QuantumOpticsBase._tp_matmul_first!(result_data, a, b_data, alpha * a.factor, beta)
    elseif a.loc == ndims(b_data)
        QuantumOpticsBase._tp_matmul_last!(result_data, a, b_data, alpha * a.factor, beta)
    else
        QuantumOpticsBase._tp_matmul_mid!(result_data, a, a.loc, b_data, alpha * a.factor, beta)
    end

    return result
end


function interaction_cpu_fallback_last!(result, op1, op2, b, α, β)
    result_parent = parent_of(result)
    b_parent = parent_of(b)

    result_host = Array(result_parent)
    b_host = Array(b_parent)

    @inbounds for row in axes(result_host, 1)
        WaveguideQED.waveguide_interaction_mul!(
            view(result_host, row, :),
            op1,
            op2,
            view(b_host, row, :),
            α,
            β,
        )
    end

    copyto!(result_parent, result_host)
    return result
end


function interaction_cpu_fallback_first!(result, op1, op2, b, α, β)
    result_parent = parent_of(result)
    b_parent = parent_of(b)

    result_host = Array(result_parent)
    b_host = Array(b_parent)

    @inbounds for col in axes(result_host, 2)
        WaveguideQED.waveguide_interaction_mul!(
            view(result_host, :, col),
            op1,
            op2,
            view(b_host, :, col),
            α,
            β,
        )
    end

    copyto!(result_parent, result_host)
    return result
end


interaction_cpu_fallback_last!(result, a::WaveguideInteraction, b, α, β) =
    interaction_cpu_fallback_last!(result, a.op1, a.op2, b, α, β)


interaction_cpu_fallback_first!(result, a::WaveguideInteraction, b, α, β) =
    interaction_cpu_fallback_first!(result, a.op1, a.op2, b, α, β)


function interaction_tmp_like(x)
    xp = parent_of(x)
    tmp_parent = similar(xp)
    fill!(tmp_parent, zero(eltype(tmp_parent)))
    return reshape(tmp_parent, size(x))
end


@inline interaction_single_index(nsteps, idx, timeindex) = 1 + (idx - 1) * nsteps + timeindex


@inline function interaction_same_offset(Nw, nsteps, idx)
    tri = (nsteps * (nsteps + 1)) ÷ 2
    return 1 + Nw * nsteps + (idx - 1) * tri
end


@inline function interaction_cross_offset(Nw, nsteps, idx1, idx2)
    tri = (nsteps * (nsteps + 1)) ÷ 2
    i = min(idx1, idx2)
    j = max(idx1, idx2)
    pair_index = (i - 1) * Nw + j - (i * (i + 1)) ÷ 2
    return 1 + Nw * nsteps + Nw * tri + (pair_index - 1) * nsteps^2
end


@inline function interaction_same_slice_raw_index(Nw, nsteps, idx, fixed_time, x)
    offset = interaction_same_offset(Nw, nsteps, idx)
    if x < fixed_time
        return offset + twophoton_index(x, nsteps, fixed_time)
    else
        return offset + twophoton_index(fixed_time, nsteps, x)
    end
end


@inline function interaction_cross_slice_raw_index(Nw, nsteps, fixed_idx, other_idx, fixed_time, x)
    offset = interaction_cross_offset(Nw, nsteps, fixed_idx, other_idx)
    if fixed_idx < other_idx
        return offset + (x - 1) * nsteps + fixed_time
    else
        return offset + (fixed_time - 1) * nsteps + x
    end
end


@inline function interaction_same_slice_match(raw_index, Nw, nsteps, idx, fixed_time)
    offset = interaction_same_offset(Nw, nsteps, idx)
    rel = raw_index - offset
    diag = twophoton_index(fixed_time, nsteps, fixed_time)

    if rel == diag
        return true, fixed_time
    elseif diag < rel <= twophoton_index(fixed_time, nsteps, nsteps)
        return true, raw_index - offset - twophoton_index(fixed_time, nsteps, 0)
    else
        is_match, x = reverse_idx(rel, nsteps, fixed_time)
        return is_match, is_match ? x : 0
    end
end


@inline function interaction_cross_slice_match(raw_index, Nw, nsteps, fixed_idx, other_idx, fixed_time)
    offset = interaction_cross_offset(Nw, nsteps, fixed_idx, other_idx)
    rel = raw_index - offset

    if fixed_idx < other_idx
        is_match, x = reverse_idx_waveguide(rel, nsteps, fixed_time)
        return is_match, is_match ? x : 0
    else
        lower = (fixed_time - 1) * nsteps + 1
        upper = lower + nsteps - 1
        is_match = lower <= rel <= upper
        return is_match, is_match ? rel - (fixed_time - 1) * nsteps : 0
    end
end


function apply_last_op_gpu!(
    result,
    op1::WaveguideCreate{B1,B2,2,idx1},
    op2::WaveguideDestroy{B1,B2,2,idx2},
    b,
    α,
    β,
) where {B1,B2,idx1,idx2}
    N, M = size(result, 1), size(result, 2)
    nsteps = op1.basis_l.nsteps
    timeindex_1 = (op1.timeindex + op1.delay - 1) % nsteps + 1
    timeindex_2 = (op2.timeindex + op2.delay - 1) % nsteps + 1
    Nw = get_number_of_waveguides(op1.basis_l)
    factor = op1.factor * op2.factor

    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_destroy_2photon_last!(
        result,
        b,
        α,
        β,
        N,
        M,
        nsteps,
        timeindex_1,
        timeindex_2,
        factor,
        Nw,
        idx1,
        idx2,
    )
    return result
end


function apply_last_op_gpu!(
    result,
    op1::WaveguideCreate{B1,B2,2,idx1},
    op2::WaveguideCreate{B1,B2,2,idx2},
    b,
    α,
    β,
) where {B1,B2,idx1,idx2}
    N, M = size(result, 1), size(result, 2)
    nsteps = op1.basis_l.nsteps
    timeindex_1 = (op1.timeindex + op1.delay - 1) % nsteps + 1
    timeindex_2 = (op2.timeindex + op2.delay - 1) % nsteps + 1
    Nw = get_number_of_waveguides(op1.basis_l)
    factor = op1.factor * op2.factor

    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_create_create_2photon_last!(
        result,
        b,
        α,
        β,
        N,
        M,
        nsteps,
        timeindex_1,
        timeindex_2,
        factor,
        Nw,
        idx1,
        idx2,
    )
    return result
end


function apply_last_op_gpu!(
    result,
    op1::WaveguideDestroy{B1,B2,2,idx1},
    op2::WaveguideCreate{B1,B2,2,idx2},
    b,
    α,
    β,
) where {B1,B2,idx1,idx2}
    N, M = size(result, 1), size(result, 2)
    nsteps = op1.basis_l.nsteps
    timeindex_1 = (op1.timeindex + op1.delay - 1) % nsteps + 1
    timeindex_2 = (op2.timeindex + op2.delay - 1) % nsteps + 1
    Nw = get_number_of_waveguides(op1.basis_l)
    factor = op1.factor * op2.factor

    blockx, blocky = 8, 32
    gridx = cld(N, blockx)
    gridy = cld(M, blocky)

    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_interaction_destroy_create_2photon_last!(
        result,
        b,
        α,
        β,
        N,
        M,
        nsteps,
        timeindex_1,
        timeindex_2,
        factor,
        Nw,
        idx1,
        idx2,
    )
    return result
end


function gpu_interaction_create_destroy_2photon_last!(
    result,
    b,
    α,
    β,
    N,
    M,
    nsteps,
    timeindex_1,
    timeindex_2,
    factor,
    Nw,
    idx1,
    idx2,
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    if !isone(β)
        @inbounds result[row, col] *= β
    end

    αfac = α * factor
    sqrt2 = sqrt(2.0)

    if col == interaction_single_index(nsteps, idx1, timeindex_1)
        @inbounds result[row, col] += αfac * b[row, interaction_single_index(nsteps, idx2, timeindex_2)]
    end

    if idx1 == idx2
        for x in 1:nsteps
            outcol = interaction_same_slice_raw_index(Nw, nsteps, idx1, timeindex_1, x)
            if col == outcol
                incol = interaction_same_slice_raw_index(Nw, nsteps, idx1, timeindex_2, x)
                in_scale = x == timeindex_2 ? sqrt2 : 1.0
                out_scale = x == timeindex_1 ? sqrt2 : 1.0
                @inbounds result[row, col] += out_scale * αfac * in_scale * b[row, incol]
            end
        end

        for k in 1:Nw
            if k == idx1
                continue
            end
            for x in 1:nsteps
                outcol = interaction_cross_slice_raw_index(Nw, nsteps, idx1, k, timeindex_1, x)
                if col == outcol
                    incol = interaction_cross_slice_raw_index(Nw, nsteps, idx1, k, timeindex_2, x)
                    @inbounds result[row, col] += αfac * b[row, incol]
                end
            end
        end
    else
        for x in 1:nsteps
            outcol_cross = interaction_cross_slice_raw_index(Nw, nsteps, idx1, idx2, timeindex_1, x)
            if col == outcol_cross
                incol_same = interaction_same_slice_raw_index(Nw, nsteps, idx2, timeindex_2, x)
                in_scale = x == timeindex_2 ? sqrt2 : 1.0
                @inbounds result[row, col] += αfac * in_scale * b[row, incol_same]
            end

            outcol_same = interaction_same_slice_raw_index(Nw, nsteps, idx1, timeindex_1, x)
            if col == outcol_same
                incol_cross = interaction_cross_slice_raw_index(Nw, nsteps, idx1, idx2, timeindex_2, x)
                out_scale = x == timeindex_1 ? sqrt2 : 1.0
                @inbounds result[row, col] += out_scale * αfac * b[row, incol_cross]
            end
        end

        for k in 1:Nw
            if k == idx1 || k == idx2
                continue
            end
            for x in 1:nsteps
                outcol = interaction_cross_slice_raw_index(Nw, nsteps, idx1, k, timeindex_1, x)
                if col == outcol
                    incol = interaction_cross_slice_raw_index(Nw, nsteps, idx2, k, timeindex_2, x)
                    @inbounds result[row, col] += αfac * b[row, incol]
                end
            end
        end
    end

    return
end


function gpu_interaction_create_create_2photon_last!(
    result,
    b,
    α,
    β,
    N,
    M,
    nsteps,
    timeindex_1,
    timeindex_2,
    factor,
    Nw,
    idx1,
    idx2,
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    if !isone(β)
        @inbounds result[row, col] *= β
    end

    αfac = α * factor
    sqrt2 = sqrt(2.0)

    if idx1 == idx2
        t1 = min(timeindex_1, timeindex_2)
        t2 = max(timeindex_1, timeindex_2)
        target = interaction_same_offset(Nw, nsteps, idx1) + twophoton_index(t1, nsteps, t2)
        if col == target
            diag_scale = timeindex_1 == timeindex_2 ? sqrt2 : 1.0
            @inbounds result[row, col] += diag_scale * αfac * b[row, 1]
        end
    else
        t1 = idx1 >= idx2 ? timeindex_1 : timeindex_2
        t2 = idx1 >= idx2 ? timeindex_2 : timeindex_1
        target = interaction_cross_offset(Nw, nsteps, idx1, idx2) + (t1 - 1) * nsteps + t2
        if col == target
            @inbounds result[row, col] += αfac * b[row, 1]
        end
    end

    return
end


function gpu_interaction_destroy_create_2photon_last!(
    result,
    b,
    α,
    β,
    N,
    M,
    nsteps,
    timeindex_1,
    timeindex_2,
    factor,
    Nw,
    idx1,
    idx2,
)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if row > N || col > M
        return
    end

    if !isone(β)
        @inbounds result[row, col] *= β
    end

    αfac = α * factor
    sqrt2 = sqrt(2.0)
    if idx1 == idx2 && timeindex_1 == timeindex_2
        if col == 1
            @inbounds result[row, col] += αfac * b[row, 1]
        end

        for j in 1:Nw
            for x in 1:nsteps
                target = interaction_single_index(nsteps, j, x)
                if col == target
                    same_site_factor = (j == idx1 && x == timeindex_1) ? 2.0 : 1.0
                    @inbounds result[row, col] += same_site_factor * b[row, target]
                end
            end
        end
    elseif col == interaction_single_index(nsteps, idx2, timeindex_2)
        @inbounds result[row, col] += αfac * b[row, interaction_single_index(nsteps, idx1, timeindex_1)]
    end

    return
end


# Apply op2 first, then op1, using the existing single-operator GPU kernels.
function interaction_apply_last_gpu!(result, op1::InteractionWaveguideOp, op2::InteractionWaveguideOp, b, α, β)
    tmp = interaction_tmp_like(b)
    apply_last_op_gpu!(tmp, op2, b, one(α), zero(β))
    apply_last_op_gpu!(result, op1, tmp, α, β)
    return result
end


function interaction_apply_first_gpu!(result, op1::InteractionWaveguideOp, op2::InteractionWaveguideOp, b, α, β)
    tmp = interaction_tmp_like(b)
    apply_first_op_gpu!(tmp, op2, b, one(α), zero(β))
    apply_first_op_gpu!(result, op1, tmp, α, β)
    return result
end


function apply_last_op_gpu!(result, op1::InteractionWaveguideOp, op2::InteractionWaveguideOp, b, α, β)
    return interaction_apply_last_gpu!(result, op1, op2, b, α, β)
end


function apply_first_op_gpu!(result, op1::InteractionWaveguideOp, op2::InteractionWaveguideOp, b, α, β)
    return interaction_apply_first_gpu!(result, op1, op2, b, α, β)
end


function apply_last_op_gpu!(result, op1, op2, b, α, β)
    return interaction_cpu_fallback_last!(result, op1, op2, b, α, β)
end


function apply_first_op_gpu!(result, op1, op2, b, α, β)
    return interaction_cpu_fallback_first!(result, op1, op2, b, α, β)
end


function QuantumOpticsBase._tp_matmul_mid!(result, a::WaveguideInteraction, loc::Integer, b, α::Number, β::Number)
    sz_b_1 = 1
    for i in 1:loc-1
        sz_b_1 *= size(b, i)
    end
    sz_b_3 = 1
    for i in loc+1:ndims(b)
        sz_b_3 *= size(b, i)
    end

    br = Base.ReshapedArray(b, (sz_b_1, size(b, loc), sz_b_3), ())
    result_r = Base.ReshapedArray(result, (sz_b_1, size(a, 1), sz_b_3), ())

    move_left = sz_b_1 < sz_b_3
    perm = move_left ? (2, 1, 3) : (1, 3, 2)

    br_p = QuantumOpticsBase._tp_matmul_get_tmp(eltype(br), ((size(br, i) for i in perm)...,), :_tp_matmul_mid_in, br)
    @strided permutedims!(br_p, br, perm)

    result_r_p = QuantumOpticsBase._tp_matmul_get_tmp(eltype(result_r), ((size(result_r, i) for i in perm)...,), :_tp_matmul_mid_out, result_r)
    @strided permutedims!(result_r_p, result_r, perm)

    if move_left
        QuantumOpticsBase._tp_matmul_first!(result_r_p, a, br_p, α, β)
    else
        QuantumOpticsBase._tp_matmul_last!(result_r_p, a, br_p, α, β)
    end

    @strided permutedims!(result_r, result_r_p, perm)
    return result
end


function QuantumOpticsBase._tp_matmul_first!(
    result::CuReshapedOrCuArray{T,N1},
    a::WaveguideInteraction,
    b::CuReshapedOrCuArray{T,N2},
    α::Number,
    β::Number,
) where {T,N1,N2}
    d_first = size(b, 1)
    d_rest = length(b) ÷ d_first
    bp = parent_of(b)
    rp = parent_of(result)

    br = reshape(bp, (d_first, d_rest))
    result_r = reshape(rp, (size(a, 1), d_rest))
    apply_first_op_gpu!(result_r, a.op1, a.op2, br, α, β)
    return result
end


function QuantumOpticsBase._tp_matmul_last!(
    result::CuReshapedOrCuArray{T,N1},
    a::WaveguideInteraction,
    b::CuReshapedOrCuArray{T,N2},
    α::Number,
    β::Number,
) where {T,N1,N2}
    d_last = size(b, ndims(b))
    d_rest = length(b) ÷ d_last
    bp = parent_of(b)
    rp = parent_of(result)

    br = reshape(bp, (d_rest, d_last))
    result_r = reshape(rp, (d_rest, d_last))
    apply_last_op_gpu!(result_r, a.op1, a.op2, br, α, β)
    return result
end

