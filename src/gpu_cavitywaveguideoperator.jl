

function mul!(result::Ket{B,A},a::CavityWaveguideAbsorption,b::Ket{B,A},alpha::Number,beta::Number) where {B<:Basis,A<:CuArray}
    # == 0) Gather dimension info
    i,j = a.loc[1], a.loc[2]
    @views b_data = reshape(b.data, QuantumOpticsBase._comp_size(basis(b)))
    @views result_data = reshape(result.data, QuantumOpticsBase._comp_size(basis(result)))

    
    nd = ndims(b_data)
    dims_b = size(b_data)
    @assert nd == ndims(result_data) "b and result should have same ndims"
    @assert length(size(result_data)) == nd
    @assert 1 ≤ i ≤ nd && 1 ≤ j ≤ nd && i != j


    if nd ==2 && i == 2
        for i in 2:size(b_data,1)
            apply_first_op_gpu!(view(result_data,i,:,:),a.op,view(b_data,i-1,:,:),sqrt(i-1)*alpha*a.factor,beta)
        end
        #rmul!(view(result.data, 1 : dims_b[1] : 1 + (dims_b[2]-1)*dims_b[1]),beta)
        rmul!(view(result_data, 1,:,:),beta)
        return result
    end
    if nd ==2 && i == 1
        for i in 2:size(b_data,2)
            apply_first_op_gpu!(view(result_data,:,i,:),a.op,view(b_data,:,i-1,:),sqrt(i-1)*alpha*a.factor,beta)
        end
        #rmul!(view(result.data,1 : 1 : dims_b[1]),beta)
        rmul!(view(result_data, :,1,:),beta)
        return result
    end

    allaxes = collect(1:length(dims_b))
    leftover = [ax for ax in allaxes if ax != i && ax != j]
    neworder = vcat(leftover, j, i)

    leftover_product = prod(dims_b[ax] for ax in leftover)
    newshape = (leftover_product, dims_b[j], dims_b[i])


    b_p = QuantumOpticsBase._tp_matmul_get_tmp(
        eltype(b_data),
        (( size(b_data,k) for k in neworder)...,),
        :_tp_matmul_mid_wg_in,
        b_data
    )
    
    permutedims!(b_p, b_data, neworder)
    @views b_r_p = reshape(b_p, newshape)


    result_p = QuantumOpticsBase._tp_matmul_get_tmp(
        eltype(result_data),
        (( size(result_data,k) for k in neworder)...,),
        :_tp_matmul_mid_wg_out,
        result_data
    )
    
    permutedims!(result_p, result_data, neworder)
    @views result_r_p = reshape(result_p, newshape)

    
    # == 5) Call either _tp_matmul_first! or _tp_matmul_last!
    #   In your code, _tp_matmul_first! expects shape (d_first, d_rest) 
    #   or something similar. So adapt as needed.
    
    for i in 2:size(b_r_p,2)
        apply_last_op_gpu!(view(result_r_p,:,i,:),a.op,view(b_r_p,:,i-1,:),sqrt(i-1)*alpha*a.factor,beta)
    end
    view(result_r_p,:,1,:) .*= beta

    # == 6) Permute back
    oldorder = invperm(neworder)
    permutedims!(result_data, result_p, oldorder)
    

    # That’s it. `result` is now updated, since `result_r` is just a reshaped view.
    return result
end


function mul!(result::Ket{B,A},a::CavityWaveguideEmission,b::Ket{B,A},alpha::Number,beta::Number) where {B<:Basis,A<:CuArray}
    # == 0) Gather dimension info
    i,j = a.loc[1], a.loc[2]
    @views b_data = reshape(b.data, QuantumOpticsBase._comp_size(basis(b)))
    @views result_data = reshape(result.data, QuantumOpticsBase._comp_size(basis(result)))

    
    nd = ndims(b_data)
    dims_b = size(b_data)
    @assert nd == ndims(result_data) "b and result should have same ndims"
    @assert length(size(result_data)) == nd
    @assert 1 ≤ i ≤ nd && 1 ≤ j ≤ nd && i != j

    if nd ==2 && i == 2
        for i in 1:size(b_data,1)-1
            apply_first_op_gpu!(view(result_data,i,:,:),a.op,view(b_data,i+1,:,:),sqrt(i)*alpha*a.factor,beta)
        end
        #rmul!(view(result.data,dims_b[1] : dims_b[1] : dims_b[1] + (dims_b[2]-1)*dims_b[1]),beta)
        rmul!(view(result_data,size(b_data,1),:,:),beta)
        return result
    end
    if nd ==2 && i == 1
        for i in 1:size(b_data,2)-1
            apply_first_op_gpu!(view(result_data,:,i,:),a.op,view(b_data,:,i+1,:),sqrt(i)*alpha*a.factor,beta)
        end
        rmul!(view(result_data,:,size(b_data,2),:),beta)
        #rmul!(view(result.data,(dims_b[2]-1)*dims_b[1] +1 : 1 :(dims_b[2]-1)*dims_b[1]+ dims_b[1]),beta)
        return result
    end

    allaxes = collect(1:length(dims_b))
    leftover = [ax for ax in allaxes if ax != i && ax != j]
    neworder = vcat(leftover, j, i)

    leftover_product = prod(dims_b[ax] for ax in leftover)
    newshape = (leftover_product, dims_b[j], dims_b[i])


    b_p = QuantumOpticsBase._tp_matmul_get_tmp(
        eltype(b_data),
        (( size(b_data,k) for k in neworder)...,),
        :_tp_matmul_mid_wg_in,
        b_data
    )
    permutedims!(b_p, b_data, neworder)
    b_r_p = @views reshape(b_p, newshape)


    result_p = QuantumOpticsBase._tp_matmul_get_tmp(
        eltype(result_data),
        (( size(result_data,k) for k in neworder)...,),
        :_tp_matmul_mid_wg_out,
        result_data
    )
    permutedims!(result_p, result_data, neworder)
    result_r_p = @views reshape(result_p, newshape)


    # == 5) Call either _tp_matmul_first! or _tp_matmul_last!
    #   In your code, _tp_matmul_first! expects shape (d_first, d_rest) 
    #   or something similar. So adapt as needed.
    
    for i in 1:size(b_r_p,2)-1
        apply_last_op_gpu!(view(result_r_p,:,i,:),a.op,view(b_r_p,:,i+1,:),sqrt(i)*alpha*a.factor,beta)
    end
    #rmul!(view(result_r_p,:,size(b_r_p,2),:),beta)
    view(result_r_p,:,size(b_r_p,2),:) .*= beta

    #for k in 1:dims_b[1]
    #    view(result_p, k + (dims_b[2]-1)*dims_b[1] : dims_b[1]*dims_b[2] : k + (dims_b[2]-1)*dims_b[1] + (dims_b[3]-1)*dims_b[1]*dims_b[2])
    #end
    #(k + (j-1)*m) : m*n : (k + (j-1)*m) + (o-1)*m*n
    # == 6) Permute back
    oldorder = invperm(neworder)
    permutedims!(result_data, result_p, oldorder)
    
    # That’s it. `result` is now updated, since `result_r` is just a reshaped view.
    return result
end
