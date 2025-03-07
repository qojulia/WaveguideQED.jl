

function QuantumOpticsBase.:_tp_matmul_first!(result::Base.ReshapedArray{T,N1,AR1,<:Tuple}, a::WaveguideOperator, b::Base.ReshapedArray{T,N2,AR2,<:Tuple}, α::Number, β::Number) where {
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
    apply_first_op_gpu!(result_r,a,br,α,β)
    result
end

#Same as _tp_matmul_first! But indexed in another way.
function QuantumOpticsBase.:_tp_matmul_last!(result::Base.ReshapedArray{T,N1,AR1,<:Tuple}, a::WaveguideOperator, b::Base.ReshapedArray{T,N2,AR2,<:Tuple}, α::Number, β::Number) where {
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
    apply_last_op_gpu!(result_r,a,br,α,β)
    result
end


function apply_last_op_gpu!(result,op::WaveguideDestroy{B,B,1,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor

    # decide block sizes
    blockx, blocky = 8, 32
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_onephoton_destroy_last!(
    result, b,alpha,beta,
    N, M,
    nsteps,
    timeindex,
    factor,idx
    )
    return result
end


function apply_first_op_gpu!(result,op::WaveguideDestroy{B,B,1,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor

    # decide block sizes
    blockx, blocky = 32, 8
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_onephoton_destroy_first!(
    result, b,alpha,beta,
    N, M,
    nsteps,
    timeindex,
    factor,idx
    )
    return result
end




function apply_last_op_gpu!(result,op::WaveguideCreate{B,B,1,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor

    # decide block sizes
    blockx, blocky = 8, 32
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_onephoton_create_last!(
    result, b,alpha,beta,
    N, M,
    nsteps,
    timeindex,
    factor,idx
    )
    return result
end


function apply_first_op_gpu!(result,op::WaveguideCreate{B,B,1,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor

    # decide block sizes
    blockx, blocky = 32, 8
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_onephoton_create_first!(
    result, b,alpha,beta,
    N, M,
    nsteps,
    timeindex,
    factor,idx
    )
    return result
end




function gpu_waveguide_onephoton_create_last!(result, b,alpha,beta, N, M, nsteps,timeindex, factor,waveguide_idx) 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    #sr = 
    #@cuprintln (row,i)

    # Bounds check
    if row > N || i > M
        return
    end
    #idx_result = (row-1) * M + i
    if !isone(beta)
        @inbounds result[row,i] *= beta
    end

    if i==1
        @inbounds result[row,timeindex+(waveguide_idx-1)*nsteps+1] += alpha*factor*b[row,i]
    end
    
    return
end


function gpu_waveguide_onephoton_create_first!(result, b,alpha,beta, N, M, nsteps,timeindex, factor,waveguide_idx) 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    # Bounds check
    if row > N || i > M
        return
    end
    #idx_result = (row-1) * M + i
    if !isone(beta)
        @inbounds result[row,i] *= beta
    end

    if row==1
        @inbounds result[timeindex+(waveguide_idx-1)*nsteps+1,i] += alpha*factor*b[row,i]
    end
    
    return
end


function gpu_waveguide_onephoton_destroy_last!(result, b,alpha,beta, N, M, nsteps,timeindex, factor,waveguide_idx) 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    #sr = 
    #@cuprintln (row,i)

    # Bounds check
    if row > N || i > M
        return
    end
    #idx_result = (row-1) * M + i
    if !isone(beta)
        @inbounds result[row,i] *= beta
    end

    if i==1
        @inbounds result[row,i] += alpha*factor*b[row,timeindex+(waveguide_idx-1)*nsteps+1]
    end
    
    return
end


function gpu_waveguide_onephoton_destroy_first!(result, b,alpha,beta, N, M, nsteps,timeindex, factor,waveguide_idx) 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    #sr = 
    #@cuprintln (row,i)

    # Bounds check
    if row > N || i > M
        return
    end
    #idx_result = (row-1) * M + i
    if !isone(beta)
        @inbounds result[row,i] *= beta
    end

    if row==1
        @inbounds result[row,i] += alpha*factor*b[timeindex+(waveguide_idx-1)*nsteps+1,i]
    end
    
    return
end



#TWOphoton part

function apply_last_op_gpu!(result,op::WaveguideDestroy{B,B,2,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor
    Nw        = get_number_of_waveguides(op.basis_l)

    # decide block sizes
    blockx, blocky = 8, 32
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_twophoton_destroy_last!(
    result, b, alpha, beta,
        N, M,
        nsteps,
        timeindex,
        factor,
        idx,     # which waveguide
        Nw       # total waveguides
    )
    return result
end


function apply_first_op_gpu!(result,op::WaveguideDestroy{B,B,2,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor
    Nw        = get_number_of_waveguides(op.basis_l)

    # decide block sizes
    blockx, blocky = 32, 8
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_twophoton_destroy_first!(
    result, b, alpha, beta,
        N, M,
        nsteps,
        timeindex,
        factor,
        idx,     # which waveguide
        Nw       # total waveguides
    )
    return result
end



function apply_last_op_gpu!(result,op::WaveguideCreate{B,B,2,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor
    Nw        = get_number_of_waveguides(op.basis_l)

    # decide block sizes
    blockx, blocky = 8, 32
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_twophoton_create_last!(
    result, b, alpha, beta,
        N, M,
        nsteps,
        timeindex,
        factor,
        idx,     # which waveguide
        Nw       # total waveguides
    )
    return result
end


function apply_first_op_gpu!(result,op::WaveguideCreate{B,B,2,idx},b,alpha, beta) where{B,idx}

    N, M = size(result, 1), size(result, 2)
    
    # read waveguide info
    nsteps    = op.basis_l.nsteps
    timeindex = (op.timeindex + op.delay - 1) % nsteps + 1
    factor    = op.factor
    Nw        = get_number_of_waveguides(op.basis_l)

    # decide block sizes
    blockx, blocky = 32, 8
    gridx  = cld(N, blockx)
    gridy  = cld(M, blocky)
    
    # Launch the kernel
    @cuda threads=(blockx, blocky) blocks=(gridx, gridy) gpu_waveguide_twophoton_create_first!(
    result, b, alpha, beta,
        N, M,
        nsteps,
        timeindex,
        factor,
        idx,     # which waveguide
        Nw       # total waveguides
    )
    return result
end


function gpu_waveguide_twophoton_destroy_last!(result, b,α, β,N, M,nsteps,timeindex,factor,waveguide_idx,Nw)
    # 2D thread indexing, 1-based or 0-based. 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    # boundary check
    if row > N || col > M
        return
    end

    # (A) scale by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # (B) single-photon part 
    if col == 1
        @inbounds result[row, col] += α*factor * b[row, timeindex + (waveguide_idx-1)*nsteps + 1]
    end

    # (C) 2-photon “same waveguide” part
    let startcol = 2 + (waveguide_idx-1)*nsteps
        endcol   = (waveguide_idx)*nsteps + 1

        if startcol <= col <= endcol
            colprime = col - startcol + 1
            offset   = Nw*nsteps + 1 + (waveguide_idx-1)*(nsteps*(nsteps+1)>>>1)

            if colprime < timeindex
                # result[row,col] += α*factor * b[row, offset + twophoton_index(colprime, timeindex)]
                result[row, col] += α*factor * b[row, offset + twophoton_index(colprime, nsteps, timeindex)]
            elseif colprime > timeindex
                result[row, col] += α*factor * b[row, offset + twophoton_index(timeindex, nsteps, colprime)]
            else
                # colprime == timeindex
                result[row, col] += sqrt(2)*α*factor * b[row, offset + twophoton_index(timeindex, nsteps, timeindex)]
            end
        end
    end

    # (D) cross-waveguide two-photon part
    
    for k in 1:Nw
        if k == waveguide_idx
            continue
        end
        startcol_k = 2 + (k-1)*nsteps
        endcol_k   = k*nsteps + 1

        if startcol_k <= col <= endcol_k
            colprime = col - startcol_k + 1
        
            i,j = min(k,waveguide_idx),max(k,waveguide_idx)
            index = (i-1)*Nw + j - (i*(i+1))÷2
            
            offset_k = 1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2

            order = (i == waveguide_idx)
            
            if order
                result[row, col] += α*factor * b[row, offset_k + (colprime-1)*nsteps + timeindex]
            else
                result[row, col] += α*factor * b[row, offset_k + (timeindex-1)*nsteps + colprime]
            end
        end
    
    end
    
    return
end






function gpu_waveguide_twophoton_destroy_first!(result, b,α, β,N, M,nsteps,timeindex,factor,waveguide_idx,Nw)
    # 2D thread indexing, 1-based or 0-based. 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    # boundary check
    if row > N || col > M
        return
    end

    # (A) scale by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end

    # (B) single-photon part 
    if row == 1
        @inbounds result[row, col] += α*factor * b[timeindex + (waveguide_idx-1)*nsteps + 1,col]
    end

    # (C) 2-photon “same waveguide” part
    let startcol = 2 + (waveguide_idx-1)*nsteps
        endcol   = (waveguide_idx)*nsteps + 1

        if startcol <= row <= endcol
            colprime = row - startcol + 1
            offset   = Nw*nsteps + 1 + (waveguide_idx-1)*(nsteps*(nsteps+1)>>>1)

            if colprime < timeindex
                # result[row,col] += α*factor * b[row, offset + twophoton_index(colprime, timeindex)]
                result[row, col] += α*factor * b[offset + twophoton_index(colprime, nsteps, timeindex),col]
            elseif colprime > timeindex
                result[row, col] += α*factor * b[offset + twophoton_index(timeindex, nsteps, colprime),col]
            else
                # colprime == timeindex
                result[row, col] += sqrt(2)*α*factor * b[offset + twophoton_index(timeindex, nsteps, timeindex),col]
            end
        end
    end

    # (D) cross-waveguide two-photon part
    
    for k in 1:Nw
        if k == waveguide_idx
            continue
        end
        startcol_k = 2 + (k-1)*nsteps
        endcol_k   = k*nsteps + 1

        if startcol_k <= row <= endcol_k
            colprime = row - startcol_k + 1
        
            i,j = min(k,waveguide_idx),max(k,waveguide_idx)
            index = (i-1)*Nw + j - (i*(i+1))÷2
            
            offset_k = 1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2

            order = (i == waveguide_idx)
            
            if order
                result[row, col] += α*factor * b[offset_k + (colprime-1)*nsteps + timeindex,col]
            else
                result[row, col] += α*factor * b[offset_k + (timeindex-1)*nsteps + colprime,col]
            end
        end
    
    end
    
    return
end


@inline function reverse_idx(C::Int,nsteps::Int,timeindex::Int)
    out = nsteps+1/2-sqrt(4*nsteps^2 - 8*C-4*nsteps + 8*timeindex + 1)
    is_int = isinteger(out)
    is_int,is_int ? Int(out) : out
end

@inline function reverse_idx_waveguide(C::Int,nsteps::Int,timeindex::Int)
    out = (C-timeindex)/nsteps+1
    is_int = isinteger(out) && (timeindex <= C <= (nsteps-1)*nsteps + timeindex)
    is_int,is_int ? Int(out) : out
end



function gpu_waveguide_twophoton_create_last!(result, b,α, β,N, M,nsteps,timeindex,factor,waveguide_idx,Nw)
    # 2D thread indexing, 1-based or 0-based. 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    # boundary check
    if row > N || col > M
        return
    end

    # (A) scale by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end
    
    # (C) 2-photon “same waveguide” part
    startcol = 2 + (waveguide_idx-1)*nsteps
    #endcol   = (waveguide_idx)*nsteps + 1
    offset   = Nw*nsteps + 1 + (waveguide_idx-1)*(nsteps*(nsteps+1)>>>1)
    
    is_int,r_idx = reverse_idx(col-offset,nsteps,timeindex)

    # (B) single-photon part 
    if col == timeindex + (waveguide_idx-1)*nsteps + 1
        @inbounds result[row, col] += α*factor * b[row, 1]
    
    #twophoton_part
    elseif col-offset == twophoton_index(timeindex, nsteps, timeindex)
        colprime = col -offset - twophoton_index(timeindex, nsteps, 0)
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
        
        for k in 1:Nw
            if k == waveguide_idx
                continue
            end
            
            
            startcol_k = 2 + (k-1)*nsteps
            
            i,j = min(k,waveguide_idx),max(k,waveguide_idx)
            index = (i-1)*Nw + j - (i*(i+1))÷2
            offset_k = 1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2
            order = (i == waveguide_idx)
            #col_idx = (order ? offset_k + (colprime-1)*nsteps + timeindex : offset_k + (timeindex-1)*nsteps + colprime)

            if order
                is_in,r_idx = reverse_idx_waveguide(col-offset_k,nsteps,timeindex)
                if is_in
                    colb = r_idx + startcol_k - 1
                    @inbounds result[row, col] += α*factor * b[row, colb]
                end
            elseif (timeindex-1)*nsteps +1 <= col-offset_k <= (timeindex-1)*nsteps + nsteps
                colb = col-offset_k - (timeindex-1)*nsteps + startcol_k - 1
                @inbounds result[row, col] += α*factor * b[row, colb]
            end
        
        end
        
    end

    return
end




function gpu_waveguide_twophoton_create_first!(result, b,α, β,N, M,nsteps,timeindex,factor,waveguide_idx,Nw)
    # 2D thread indexing, 1-based or 0-based. 
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    col   = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    # boundary check
    if row > N || col > M
        return
    end

    # (A) scale by β
    if !isone(β)
        @inbounds result[row, col] *= β
    end
    
    # (C) 2-photon “same waveguide” part
    startcol = 2 + (waveguide_idx-1)*nsteps
    #endcol   = (waveguide_idx)*nsteps + 1
    offset   = Nw*nsteps + 1 + (waveguide_idx-1)*(nsteps*(nsteps+1)>>>1)
    
    is_int,r_idx = reverse_idx(row-offset,nsteps,timeindex)

    # (B) single-photon part 
    if row == timeindex + (waveguide_idx-1)*nsteps + 1
        @inbounds result[row, col] += α*factor * b[1, col]
    
    #twophoton_part
    elseif row-offset == twophoton_index(timeindex, nsteps, timeindex)
        colprime = row -offset - twophoton_index(timeindex, nsteps, 0)
        colb = colprime + startcol - 1
        @inbounds result[row, col] += sqrt(2)*α*factor * b[colb, col]

    elseif twophoton_index(timeindex, nsteps, timeindex)+1 <= row-offset <= twophoton_index(timeindex, nsteps, nsteps)
        colprime = row -offset - twophoton_index(timeindex, nsteps, 0)
        colb = colprime + startcol - 1
        
        @inbounds result[row, col] += α*factor * b[colb, col]

    elseif is_int
        colb = r_idx + startcol - 1
        @inbounds result[row, col] += α*factor * b[colb, col]
    else
        
        for k in 1:Nw
            if k == waveguide_idx
                continue
            end
            
            
            startcol_k = 2 + (k-1)*nsteps
            
            i,j = min(k,waveguide_idx),max(k,waveguide_idx)
            index = (i-1)*Nw + j - (i*(i+1))÷2
            offset_k = 1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2
            order = (i == waveguide_idx)
            #col_idx = (order ? offset_k + (colprime-1)*nsteps + timeindex : offset_k + (timeindex-1)*nsteps + colprime)

            if order
                is_in,r_idx = reverse_idx_waveguide(row-offset_k,nsteps,timeindex)
                if is_in
                    colb = r_idx + startcol_k - 1
                    @inbounds result[row, col] += α*factor * b[colb, col]
                end
            elseif (timeindex-1)*nsteps +1 <= row-offset_k <= (timeindex-1)*nsteps + nsteps
                colb = row-offset_k - (timeindex-1)*nsteps + startcol_k - 1
                @inbounds result[row, col] += α*factor * b[colb, col]
            end
        
        end
        
    end

    return
end
