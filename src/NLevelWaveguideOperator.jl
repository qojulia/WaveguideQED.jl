#abstract type  NWaveguideOperator{BL,BR} <: AbstractOperator{BL,BR} end


"""
    NLevelWaveguideOperator{B1,B2} <: NWaveguideOperator{BL,BR}

Structure for fast simultaneous NLevelTransition and Waveguide operator
"""
mutable struct NLevelWaveguideOperator{B1,B2} <: AbstractOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    op::AbstractOperator
    loc
    n_to::Int
    n_from::Int
    indexing::WaveguideIndexing
    function NLevelWaveguideOperator(basis_l::B1,basis_r::B2,factor::ComplexF64,op::AbstractOperator,loc,n_to,n_from) where{B1,B2}
        new{B1,B2}(basis_l,basis_r,factor,op,loc,n_to,n_from,WaveguideIndexing(basis_l,loc))
    end
end

@inline function set_time!(o::NLevelWaveguideOperator, t::Number)
    set_time!(o.op,t)
end

function Base.:eltype(x::NLevelWaveguideOperator) typeof(x.factor) end
function Base.:copy(x::NLevelWaveguideOperator)
    NLevelWaveguideOperator(x.basis_l,x.basis_r,x.factor,x.op,x.loc,x.n_to,x.n_from)
end

Base.:*(x::NLevelWaveguideOperator{BL,BR},y::NLevelWaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor*y.factor)
Base.:*(x::NLevelWaveguideOperator{BL,BR},y::WaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor*y.factor)
Base.:*(x::WaveguideOperator{BL,BR},y::NLevelWaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor*y.factor)

Base.:*(x::Operator{BL,BR},y::NLevelWaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),y.factor)
Base.:*(x::NLevelWaveguideOperator{BL,BR},y::Operator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor)

#Method for multiplying, which updates factor in the operator.
function Base.:*(a::Number,b::NLevelWaveguideOperator)
    out = copy(b)
    out.factor=out.factor*a
    out
end
Base.:*(b::NLevelWaveguideOperator,a::Number)=*(a,b)


function get_nlevel_operator(a::NLevelWaveguideOperator)
    transition(a.basis_l.bases[a.loc[2]],a.n_to,a.n_from)
end



_is_transition(a::Operator,i,j) = _is_transition(a.data,basis(a),i,j)
function _is_transition(a::AbstractArray,b::Basis,i,j) 
    # Get the transition operator
    trans_op = QuantumOptics.transition(b,i,j).data
    
    # Find non-zero elements in both operators
    non_zero_a = findall(!iszero, a)
    non_zero_trans = findall(!iszero, trans_op)
    
    # If the non-zero patterns don't match, it's not a transition operator
    if non_zero_a != non_zero_trans
        return false,1.0
    end
    
    # If there are no non-zero elements, both are zero operators
    if isempty(non_zero_a)
        return true,1.0
    end
    
    # Get the scaling factor from the first non-zero element
    scale_factor = a[first(non_zero_a)] / trans_op[first(non_zero_a)]
    
    # Check if all non-zero elements are scaled by the same factor
    return all(idx -> isapprox(a[idx], trans_op[idx] * scale_factor), non_zero_a),scale_factor
end

function QuantumOpticsBase.:+(a::NLevelWaveguideOperator,b::NLevelWaveguideOperator)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum(a) +LazySum(b)
end
function QuantumOpticsBase.:-(a::NLevelWaveguideOperator,b::NLevelWaveguideOperator)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum(a) - LazySum(b) 
end
function QuantumOpticsBase.:-(a::NLevelWaveguideOperator)
    out = copy(a)
    out.factor = -a.factor
    return out
end



"""
    identityoperator(a::NLevelWaveguideOperator)

Return identityoperator(a.basis_l).
"""
function QuantumOpticsBase.identityoperator(a::NLevelWaveguideOperator)
    identityoperator(a.basis_l)
end
function QuantumOpticsBase.identityoperator(::Type{T}, b1::Basis, b2::Basis) where {T<:NLevelWaveguideOperator}
    @assert b1==b2
    identityoperator(b1)
end


"""
    tensor(a::AbstractOperator,b::CavityWaveguideAbsorption)
    tensor(a::CavityWaveguideAbsorption,b::AbstractOperator)
    tensor(a::AbstractOperator,b::CavityWaveguideEmission)
    tensor(a::CavityWaveguideEmission,b::AbstractOperator)

Methods for tensorproducts between QuantumOptics.jl operators and [`NLevelWaveguideOperator`](@ref).
"""
function tensor(a::AbstractOperator,b::NLevelWaveguideOperator)
    btotal = tensor(a.basis_l,b.basis_r)
    if _is_identity(a)
        NLevelWaveguideOperator(btotal,btotal,b.factor,b.op,b.loc .+length(basis(a).shape),b.n_to,b.n_from)
    else
        sorted_idx = sortperm([1,b.loc[1]+1,b.loc[2]+1])
        LazyTensor(btotal,btotal,[1,b.loc[1]+1,b.loc[2]+1][sorted_idx],(a,b.op,get_nlevel_operator(b))[sorted_idx])
    end
end
function tensor(a::NLevelWaveguideOperator,b::AbstractOperator)
    btotal = tensor(a.basis_l,b.basis_r)
    if _is_identity(b)
        NLevelWaveguideOperator(btotal,btotal,a.factor,a.op,a.loc,a.n_to,a.n_from)
    else
        sorted_idx = sortperm([a.loc[1]+1,a.loc[2]+1,length(btotal.shape)])
        LazyTensor(btotal,btotal,[a.loc[1]+1,a.loc[2]+1,length(btotal.shape)][sorted_idx],(a.op,get_nlevel_operator(a),b)[sorted_idx])
    end
end
function tensor(a::T,b::Operator{BL,BR,F}) where {BL<:NLevelBasis,BR<:NLevelBasis,F,T<:WaveguideOperatorT}
    if _is_identity(b)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,[a.indices...],(a.operators...,),a.factor)
    end
    N_basis = basis(b)
    for i in 1:N_basis.N
        for j in 1:N_basis.N
            truth_value, prefactor = _is_transition(b,i,j)
            if truth_value
                btotal = basis(a) ⊗ basis(b)
                return NLevelWaveguideOperator(btotal,btotal,a.factor*prefactor,a.operators[1],[a.indices[1],length(basis(a).shape)+1],i,j) 
            end
        end
    end
    a ⊗ LazyTensor(b.basis_l,b.basis_r,[1],(b,),1)
end
function tensor(a::Operator{BL,BR,F},b::T) where {BL<:NLevelBasis,BR<:NLevelBasis,F,T<:WaveguideOperatorT}
    if _is_identity(a)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,[1,b.loc[1]+1,b.loc[2]+1][sorted_idx],(a,b.op,get_nlevel_operator(b))[sorted_idx])
    end
    N_basis = basis(a)
    for i in 1:N_basis.N
        for j in 1:N_basis.N
            truth_value, prefactor = _is_transition(a,i,j)
            if truth_value
                btotal = basis(a) ⊗ basis(b)
                return NLevelWaveguideOperator(btotal,btotal,b.factor*prefactor,b.operators[1],[b.indices[1],1],i,j) 
            end
        end
    end
    LazyTensor(a.basis_l,a.basis_r,[1],(a,),1) ⊗ b
end
function tensor(a::T,b::Operator{BL,BR,F})  where {BL<:NLevelBasis,BR<:NLevelBasis,F,T<:WaveguideOperator}
    if _is_identity(b)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,(1,),(a,))
    end
    N_basis = basis(b)
    for i in 1:N_basis.N
        for j in 1:N_basis.N
            truth_value, prefactor = _is_transition(a,i,j)
            if truth_value
                btotal = basis(a) ⊗ basis(b)
                return NLevelWaveguideOperator(btotal,btotal,prefactor,a,[1,2],i,j) 
            end
        end
    end
    btotal = basis(a) ⊗ basis(b)
    return LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
end
function tensor(a::Operator{BL,BR,F},b::T) where {BL<:NLevelBasis,BR<:NLevelBasis,F,T<:WaveguideOperator}
    if _is_identity(a)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,(length(basis(a).shape)+1,),(b,))
    end
    N_basis = basis(a)
    for i in 1:N_basis.N
        for j in 1:N_basis.N
            truth_value, prefactor = _is_transition(a,i,j)
            if truth_value
                btotal = basis(a) ⊗ basis(b)
                return NLevelWaveguideOperator(btotal,btotal,prefactor,b,[2,1],i,j) 
            end
        end
    end
    btotal = basis(a) ⊗ basis(b)
    LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
end


#TO DO: CLEAN UP
#Currently indexing is very complicated using the CavityIndexing structure and could possible be done smoother and faster.
"""
    mul!(result::Ket{B1}, a::CavityWaveguideEmission, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis}
    mul!(result::Ket{B1}, a::CavityWaveguideAbsorption, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis}
    
Fast in-place multiplication of operators/state vectors. Updates `result` as `result = alpha*a*b + beta*result`. `a` is a [`NLevelWaveguideOperator`](@ref).
"""
function mul!(result::Ket{B1,A1}, a::NLevelWaveguideOperator, b::Ket{B2,A2}, alpha, beta) where {B1<:Basis,B2<:Basis, A1<:AbstractArray, A2<:AbstractArray}
    dims = basis(result).shape
    i,j = a.loc[1],a.loc[2]
    beta1 = beta
    if iszero(beta)
        fill!(result.data,beta)
        beta1 = one(beta)
    end
    if length(dims) == 2
        loop_transition_ax!(result.data,a,b.data,alpha,beta1,a.indexing,a.n_to,a.n_from)                     
    elseif length(dims) == 3
        loop_transition_third_axis!(result.data,a,b.data,alpha,beta1,a.indexing,1,a.n_to,a.n_from)     
        a.indexing.idx_vec1[j] = dims[j]
        for k in 1:dims[j]
            if i !=a.n_to || k!=a.n_from 
                loop_rmul_axis!(result.data,a,b.data,alpha,beta1,a.indexing,k)
            end
        end
    else
        iterate_over_iter!(result.data,a,b.data,alpha,beta1,a.indexing,1,loop_transition_third_axis!,a.n_to,a.n_from)
        for k in 1:dims[j]
            if i !=a.n_to || k!=a.n_from 
                a.indexing.idx_vec1[j] = k
                iterate_over_iter!(result.data,a,b.data,alpha,beta1,a.indexing,1,loop_rmul_axis!)
            end
        end 
    end
    return result
end

function loop_transition_ax!(result,a,b,alpha,beta,indexing::WaveguideIndexing,k,l)
    i,j = a.loc[1],a.loc[2]
    indexing.idx_vec1[j] = k
    indexing.idx_vec2[j] = l
    li1 = linear_index(indexing.ndims,indexing.idx_vec1,indexing.strides)
    li2 = linear_index(indexing.ndims,indexing.idx_vec2,indexing.strides)  
    waveguide_mul!(view(result,li1:indexing.strides[a.loc[1]]:li1+(indexing.end_idx-1)*indexing.strides[a.loc[1]]),a.op,view(b,li2:indexing.strides[a.loc[1]]:li2+(indexing.end_idx-1)*indexing.strides[a.loc[1]]),a.factor*alpha,beta)
end
function loop_transition_third_axis!(result,a,b,alpha,beta,indexing::WaveguideIndexing,idx,k,l)
    for j in 1:indexing.iter[idx]
        @inbounds indexing.idx_vec1[indexing.range_idx[idx]] = j
        @inbounds indexing.idx_vec2[indexing.range_idx[idx]] = j
        loop_transition_ax!(result,a,b,alpha,beta,indexing,k,l)
    end
end
function iterate_over_iter!(result,b,a,alpha,beta,indexing::WaveguideIndexing,idx,f::Function,k,m)
    if length(indexing.iter) == idx
        f(result,b,a,alpha,beta,indexing,idx,k,m)
    else
    for j in 1:indexing.iter[idx]
        @inbounds indexing.idx_vec1[indexing.range_idx[idx]] = j
        @inbounds indexing.idx_vec2[indexing.range_idx[idx]] = j
        iterate_over_iter!(result,b,a,alpha,beta,indexing,idx+1,f,k,m)
    end
    end
end



function get_waveguide_operators(op::NLevelWaveguideOperator)
    [op.op]
end
function set_waveguidetimeindex!(op::NLevelWaveguideOperator,timeindex)
    op.op.timeindex = timeindex
end
