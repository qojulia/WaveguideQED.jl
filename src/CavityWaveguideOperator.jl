"""
Abstract type used for operators on acting on a combined [`WaveguideBasis`](@ref) and cavity basis (`FockBasis`)
"""
abstract type  CavityWaveguideOperator{BL,BR} <: AbstractOperator{BL,BR} end

#TO DO: CLEAN UP
#Indexing structure used to loop over state. Needs cleanup and can possible be removed with the addition of eachslice(a,dims=(1,3,...)) in julia 1.9
struct WaveguideIndexing{N}
    ndims::Int
    k_idx::Int
    end_idx::Int
    strides::Vector
    idx_vec1::Vector
    idx_vec2::Vector
    range_idx::Vector
    iter::Tuple
    function WaveguideIndexing(end_idx,k_idx,strides::Vector,idx_vec1::Vector,idx_vec2::Vector,range_idx::Vector,iter::Tuple)
        new{length(iter)}(length(strides),k_idx,end_idx,strides,idx_vec1,idx_vec2,range_idx,iter)
    end
end
function WaveguideIndexing(b::Basis,loc)
    dims = b.shape
    alldims = [1:length(dims)...]
    i,j = loc[1],loc[2]
    exclude_dims = [i,j]
    otherdims = setdiff(alldims, exclude_dims)
    iter = (dims[otherdims]...,)

    idx_vec1=ones(Int64,length(dims))
    idx_vec2=ones(Int64,length(dims))

    range_idx = Vector{Int}(undef, length(iter))
    range_idx .= otherdims

    n = length(dims)
    prev_prod = 1
    strides = ones(Int64,n)
    for i = 2:n
        @inbounds prev_prod *= dims[i-1]
        @inbounds strides[i] = prev_prod
    end
    return WaveguideIndexing(dims[i],dims[j],strides,idx_vec1,idx_vec2,range_idx,iter)
end

"""
    CavityWaveguideAbsorption{B1,B2} <: CavityWaveguideOperator{B1,B2}

Structure for fast simultaneous Cavity creation and Waveguide annihilation operator
"""
mutable struct CavityWaveguideAbsorption{B1,B2} <: CavityWaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    op::AbstractOperator
    loc
    indexing::WaveguideIndexing
    function CavityWaveguideAbsorption(basis_l::B1,basis_r::B2,factor::ComplexF64,op::AbstractOperator,loc) where{B1,B2}
        new{B1,B2}(basis_l,basis_r,factor,op,loc,WaveguideIndexing(basis_l,loc))
    end
end

"""
    CavityWaveguideEmission{B1,B2} <: CavityWaveguideOperator{B1,B2}

Structure for fast simultaneous cavity annihilation and waveguide creation operator
"""
mutable struct CavityWaveguideEmission{B1,B2} <: CavityWaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    op::AbstractOperator
    loc
    indexing::WaveguideIndexing
    function CavityWaveguideEmission(basis_l::B1,basis_r::B2,factor::ComplexF64,op::AbstractOperator,loc) where{B1,B2}
        new{B1,B2}(basis_l,basis_r,factor,op,loc,WaveguideIndexing(basis_l,loc))
    end
end

@inline function set_time!(o::CavityWaveguideOperator, t::Number)
    set_time!(o.op,t)
end

function Base.:eltype(x::CavityWaveguideOperator) typeof(x.factor) end
function Base.:copy(x::CavityWaveguideEmission)
    CavityWaveguideEmission(x.basis_l,x.basis_r,x.factor,x.op,x.loc)
end
function Base.:copy(x::CavityWaveguideAbsorption)
    CavityWaveguideAbsorption(x.basis_l,x.basis_r,x.factor,x.op,x.loc)
end

Base.:*(x::CavityWaveguideOperator{BL,BR},y::CavityWaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor*y.factor)
Base.:*(x::CavityWaveguideOperator{BL,BR},y::WaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor*y.factor)
Base.:*(x::WaveguideOperator{BL,BR},y::CavityWaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor*y.factor)

Base.:*(x::Operator{BL,BR},y::CavityWaveguideOperator{BL,BR}) where {BL,BR} = LazyProduct((x,y),y.factor)
Base.:*(x::CavityWaveguideOperator{BL,BR},y::Operator{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor)

function Base.:*(x::Operator{BL,BR},y::WaveguideOperator{BL,BR}) where {BL,BR} 
    isa(x,QuantumOpticsBase.EyeOpType) && return y
    LazyProduct((x,y),y.factor)
end
function Base.:*(x::WaveguideOperator{BL,BR},y::Operator{BL,BR}) where {BL,BR}
    isa(y,QuantumOpticsBase.EyeOpType) && return x
    LazyProduct((x,y),x.factor)
end

const WaveguideOperatorT = LazyTensor{B,B,F,V,T} where {B,F,V,T<:Tuple{Union{WaveguideDestroy,WaveguideCreate}}}

#Method for multiplying, which updates factor in the operator.
function Base.:*(a::Number,b::CavityWaveguideOperator)
    out = copy(b)
    out.factor=out.factor*a
    out
end
Base.:*(b::CavityWaveguideOperator,a::Number)=*(a,b)


function get_cavity_operator(a::CavityWaveguideEmission)
    destroy(a.basis_l.bases[a.loc[2]])
end
function get_cavity_operator(a::CavityWaveguideAbsorption)
    create(a.basis_l.bases[a.loc[2]])
end

function QuantumOpticsBase.:+(a::CavityWaveguideOperator,b::CavityWaveguideOperator)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum(a) +LazySum(b)
end
function QuantumOpticsBase.:-(a::CavityWaveguideOperator,b::CavityWaveguideOperator)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum(a) - LazySum(b) 
end
function QuantumOpticsBase.:-(a::CavityWaveguideOperator)
    out = copy(a)
    out.factor = -a.factor
    return out
end

"""
    absorption(b1::WaveguideBasis{T},b2::FockBasis) where T
    absorption(b1::FockBasis,b2::WaveguideBasis{T}) where T

Create [`CavityWaveguideAbsorption`](@ref) that applies `create(b::FockBasis)` on `FockBasis` and destroy(b::WaveguideBasis{T}) on [`WaveguideBasis{T}`](@ref).  
"""
function absorption(b1::WaveguideBasis{T,Nw},b2::FockBasis,idx;factor=1.0) where {T,Nw}
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(factor),destroy(b1,idx),[1,2])
end
function absorption(b1::FockBasis,b2::WaveguideBasis{T,Nw},idx;factor=1.0) where {T,Nw}
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(factor),destroy(b2,idx),[2,1])
end
function absorption(b1::WaveguideBasis{T},b2::FockBasis;factor=1.0) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(factor),destroy(b1),[1,2])
end
function absorption(b1::FockBasis,b2::WaveguideBasis{T};factor=1.0) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(factor),destroy(b2),[2,1])
end
function absorption(a::T,b::Operator;factor=1.0) where T<: WaveguideOperatorT
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideAbsorption(btotal,btotal,complex(factor)*a.factor,a.operators[1],[a.indices[1],a.indices[1]+1])
end
function absorption(a::Operator,b::T;factor=1.0) where T<: WaveguideOperatorT
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideAbsorption(btotal,btotal,complex(factor)*b.factor,b.operators[1],[b.indices[1]+1,1])
end
function absorption(a::T,b::Operator;factor=1.0) where T<: WaveguideOperator
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideAbsorption(btotal,btotal,complex(factor),a,[1,2])
end
function absorption(a::Operator,b::T;factor=1.0) where T<: WaveguideOperator
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideAbsorption(btotal,btotal,complex(factor),b,[2,1])
end
function absorption(a::CompositeBasis,b::T,i::Int;factor=1.0) where T<: Union{WaveguideOperator,WaveguideOperatorT}
    btotal = tensor(a,basis(b))
    CavityWaveguideAbsorption(btotal,btotal,complex(factor),get_waveguide_operators(b)[1],[get_waveguide_location(basis(b))+length(a.shape),i])
end
function absorption(a::T,b::CompositeBasis,i::Int;factor=1.0) where T<: Union{WaveguideOperator,WaveguideOperatorT}
    btotal = tensor(basis(a),b)
    CavityWaveguideAbsorption(btotal,btotal,complex(factor),get_waveguide_operators(a)[1],[1,length(basis(a).shape)+i])
end

"""
    emission(b1::WaveguideBasis{T},b2::FockBasis) where T
    emission(b1::FockBasis,b2::WaveguideBasis{T}) where T

Create [`CavityWaveguideEmission`](@ref) that applies `destroy(b::FockBasis)` on `FockBasis` and create(b::WaveguideBasis{T}) on [`WaveguideBasis{T}`](@ref).  
"""
function emission(b1::WaveguideBasis{T,Nw},b2::FockBasis,idx;factor=1.0) where {T,Nw}
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(factor),create(b1,idx),[1,2])
end
function emission(b1::FockBasis,b2::WaveguideBasis{T,Nw},idx;factor=1.0) where {T,Nw}
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(factor),create(b2,idx),[2,1])
end
function emission(b1::WaveguideBasis{T,1},b2::FockBasis;factor=1.0) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(factor),create(b1),[1,2])
end
function emission(b1::FockBasis,b2::WaveguideBasis{T,1};factor=1.0) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(factor),create(b2),[2,1])
end
function emission(a::T,b::Operator;factor=1.0) where T<: WaveguideOperatorT
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideEmission(btotal,btotal,complex(factor)*a.factor,a.operators[1],[a.indices[1],length(basis(a).shape)+1])
end
function emission(a::Operator,b::T;factor=1.0) where T<: WaveguideOperatorT
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideEmission(btotal,btotal,complex(factor)*b.factor,b.operators[1],[b.indices[1]+1,1])
end
function emission(a::T,b::Operator;factor=1.0) where T<: WaveguideOperator
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideEmission(btotal,btotal,complex(factor),a,[1,2])
end
function emission(a::Operator,b::T;factor=1.0) where T<: WaveguideOperator
    btotal = tensor(basis(a),basis(b))
    CavityWaveguideEmission(btotal,btotal,complex(factor),b,[2,1])
end
function emission(a::CompositeBasis,b::T,i::Int;factor=1.0) where T<: Union{WaveguideOperator,WaveguideOperatorT}
    btotal = tensor(a,basis(b))
    CavityWaveguideEmission(btotal,btotal,complex(factor),get_waveguide_operators(b)[1],[get_waveguide_location(basis(b))+length(a.shape),i])
end
function emission(a::T,b::CompositeBasis,i::Int;factor=1.0) where T<: Union{WaveguideOperator,WaveguideOperatorT}
    btotal = tensor(basis(a),b)
    CavityWaveguideEmission(btotal,btotal,complex(factor),get_waveguide_operators(a)[1],[get_waveguide_location(basis(a)),length(basis(a).shape)+i])
end





"""
    identityoperator(a::CavityWaveguideOperator)

Return identityoperator(a.basis_l).
"""
function QuantumOpticsBase.identityoperator(a::CavityWaveguideOperator)
    identityoperator(a.basis_l)
end
function QuantumOpticsBase.identityoperator(::Type{T}, b1::Basis, b2::Basis) where {T<:CavityWaveguideOperator}
    @assert b1==b2
    identityoperator(b1)
end


_is_identity(a::Operator) = _is_identity(a.data,basis(a))
_is_identity(a::AbstractArray,b::Basis) = isapprox(a,identityoperator(b).data)

_is_destroy(a::Operator) = _is_destroy(a.data,basis(a))
_is_create(a::Operator) = _is_create(a.data,basis(a))
function _is_destroy(a::AbstractArray,b::FockBasis) 
    # Get the destroy operator
    destroy_op = destroy(b).data
    
    # Find non-zero elements in both operators
    non_zero_a = findall(!iszero, a)
    non_zero_destroy = findall(!iszero, destroy_op)
    
    # If the non-zero patterns don't match, it's not a destroy operator
    if non_zero_a != non_zero_destroy
        return false, 1.0
    end
    
    # If there are no non-zero elements, both are zero operators
    if isempty(non_zero_a)
        return true, 1.0
    end
    
    # Get the scaling factor from the first non-zero element
    scale_factor = a[first(non_zero_a)] / destroy_op[first(non_zero_a)]
    
    # Check if all non-zero elements are scaled by the same factor
    return all(idx -> isapprox(a[idx], destroy_op[idx] * scale_factor), non_zero_a), scale_factor
end

function _is_create(a::AbstractArray,b::FockBasis) 
    # Get the create operator
    create_op = create(b).data
    
    # Find non-zero elements in both operators
    non_zero_a = findall(!iszero, a)
    non_zero_create = findall(!iszero, create_op)
    
    # If the non-zero patterns don't match, it's not a create operator
    if non_zero_a != non_zero_create
        return false, 1.0
    end
    
    # If there are no non-zero elements, both are zero operators
    if isempty(non_zero_a)
        return true, 1.0
    end
    
    # Get the scaling factor from the first non-zero element
    scale_factor = a[first(non_zero_a)] / create_op[first(non_zero_a)]
    
    # Check if all non-zero elements are scaled by the same factor
    return all(idx -> isapprox(a[idx], create_op[idx] * scale_factor), non_zero_a), scale_factor
end




"""
    tensor(a::AbstractOperator,b::CavityWaveguideAbsorption)
    tensor(a::CavityWaveguideAbsorption,b::AbstractOperator)
    tensor(a::AbstractOperator,b::CavityWaveguideEmission)
    tensor(a::CavityWaveguideEmission,b::AbstractOperator)

Methods for tensorproducts between QuantumOptics.jl operators and [`CavityWaveguideOperator`](@ref).
"""
function tensor(a::AbstractOperator,b::CavityWaveguideAbsorption)
    btotal = tensor(a.basis_l,b.basis_r)
    if _is_identity(a)
        CavityWaveguideAbsorption(btotal,btotal,b.factor,b.op,b.loc .+length(basis(a).shape))
    else
        sorted_idx = sortperm([1,b.loc[1]+1,b.loc[2]+1])
        LazyTensor(btotal,btotal,[1,b.loc[1]+1,b.loc[2]+1][sorted_idx],(a,b.op,get_cavity_operator(b))[sorted_idx])
    end
end
function tensor(a::CavityWaveguideAbsorption,b::AbstractOperator)
    btotal = tensor(a.basis_l,b.basis_r)
    if _is_identity(b)
        CavityWaveguideAbsorption(btotal,btotal,a.factor,a.op,a.loc)
    else
        sorted_idx = sortperm([a.loc[1]+1,a.loc[2]+1,length(btotal.shape)])
        LazyTensor(btotal,btotal,[a.loc[1]+1,a.loc[2]+1,length(btotal.shape)][sorted_idx],(a.op,get_cavity_operator(a),b)[sorted_idx])
    end
end
function tensor(a::AbstractOperator,b::T) where {T<:CavityWaveguideEmission}
    btotal = tensor(a.basis_l,b.basis_r)
    if _is_identity(a)
        CavityWaveguideEmission(btotal,btotal,b.factor,b.op,b.loc .+length(basis(a).shape))
    else
        sorted_idx = sortperm([1,b.loc[1]+1,b.loc[2]+1])
        LazyTensor(btotal,btotal,[1,b.loc[1]+1,b.loc[2]+1][sorted_idx],(a,b.op,get_cavity_operator(b))[sorted_idx])
    end
end
function tensor(a::T,b::AbstractOperator) where {T<:CavityWaveguideEmission}
    btotal = tensor(a.basis_l,b.basis_r)
    if _is_identity(b)
        CavityWaveguideEmission(btotal,btotal,a.factor,a.op,a.loc)
    else
        sorted_idx = sortperm([a.loc[1]+1,a.loc[2]+1,length(btotal.shape)])
        LazyTensor(btotal,btotal,[a.loc[1]+1,a.loc[2]+1,length(btotal.shape)][sorted_idx],(a.op,get_cavity_operator(a),b)[sorted_idx])
    end
end

function tensor(a::T,b::Operator{BL,BR,F}) where {BL<:FockBasis,BR<:FockBasis,F,T<:WaveguideOperatorT}
    if _is_identity(b)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,[a.indices...],(a.operators...,),a.factor)
    end
    truth_value, prefactor = _is_destroy(b.data,basis(b))
    if truth_value
        emission(a,b;factor=prefactor)
    else
        truth_value, prefactor = _is_create(b.data,basis(b))
        if truth_value
            absorption(a,b;factor=prefactor)
        else
            a ⊗ LazyTensor(b.basis_l,b.basis_r,[1],(b,),1)
        end
    end
end
function tensor(a::Operator{BL,BR,F},b::T) where {BL<:FockBasis,BR<:FockBasis,F,T<:WaveguideOperatorT}
    if _is_identity(a)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,[b.indices...].+1 ,(b.operators...,),b.factor)
    end
    truth_value, prefactor = _is_destroy(a.data,basis(a))
    if truth_value
        emission(a,b;factor=prefactor)
    else
        truth_value, prefactor = _is_create(a.data,basis(a))
        if truth_value
            absorption(a,b;factor=prefactor)
        else
            LazyTensor(a.basis_l,a.basis_r,[1],(a,),1) ⊗ b
        end
    end
end
function tensor(a::T,b::Operator{BL,BR,F})  where {BL<:FockBasis,BR<:FockBasis,F,T<:WaveguideOperator}
    if _is_identity(b)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,(1,),(a,))
    end
    truth_value, prefactor = _is_destroy(b.data,basis(b))
    if truth_value
        emission(a,b;factor=prefactor)
    else
        truth_value, prefactor = _is_create(b.data,basis(b))
        if truth_value
            absorption(a,b;factor=prefactor)
        else
            btotal = basis(a) ⊗ basis(b)
            LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
        end
    end
end
function tensor(a::Operator{BL,BR,F},b::T) where {BL<:FockBasis,BR<:FockBasis,F,T<:WaveguideOperator}
    if _is_identity(a)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,(length(basis(a).shape)+1,),(b,))
    end
    truth_value, prefactor = _is_destroy(a.data,basis(a))
    if truth_value
        emission(a,b;factor=prefactor)
    else
        truth_value, prefactor = _is_create(a.data,basis(a))
        if truth_value
            absorption(a,b;factor=prefactor)
        else
            btotal = basis(a) ⊗ basis(b)
            LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
        end
    end
end


function _is_destroy(data,basis)
    0, 1.0
end

function _is_create(data,basis)
    0, 1.0
end

#Used to construct CavityWaveguideOperators from LazyTensors or CompositeBasis

function _is_destroy(data::AbstractArray,basis::CompositeBasis)
    N = length(basis.shape)
    ind = zeros(N)
    for k = 1:N
        ind .= 0
        ind[k] = 1
        destroy_op = tensor([ (i==0 || !isa(basis.bases[j],FockBasis)) ? identityoperator(basis.bases[j]) : destroy(basis.bases[j]) for (j,i) in enumerate(ind)]...).data
        
        # Check if this is the same destroy operator
        if isequal(data, destroy_op)
            return k, 1.0
        end
        
        # Check if it's a scaled destroy operator
        non_zero_data = findall(!iszero, data)
        non_zero_destroy = findall(!iszero, destroy_op)
        
        if non_zero_data == non_zero_destroy && !isempty(non_zero_data)
            scale_factor = data[first(non_zero_data)] / destroy_op[first(non_zero_data)]
            if all(idx -> isapprox(data[idx], destroy_op[idx] * scale_factor), non_zero_data)
                return k, scale_factor
            end
        end
    end
    return 0, 1.0
end
function _is_create(data::AbstractArray,basis::CompositeBasis)
    N = length(basis.shape)
    ind = zeros(N)
    for k = 1:N
        ind .= 0
        ind[k] = 1
        create_op = tensor([(i==0 || !isa(basis.bases[j],FockBasis)) ? identityoperator(basis.bases[j]) : create(basis.bases[j]) for (j,i) in enumerate(ind)]...).data
        
        # Check if this is the same create operator
        if isequal(data, create_op)
            return k, 1.0
        end
        
        # Check if it's a scaled create operator
        non_zero_data = findall(!iszero, data)
        non_zero_create = findall(!iszero, create_op)
        
        if non_zero_data == non_zero_create && !isempty(non_zero_data)
            scale_factor = data[first(non_zero_data)] / create_op[first(non_zero_data)]
            if all(idx -> isapprox(data[idx], create_op[idx] * scale_factor), non_zero_data)
                return k, scale_factor
            end
        end
    end
    return 0, 1.0
end

# Add this function before the tensor implementations

_is_transition(a,b) = 0, 0, 0,1.0

_is_transition(a::Operator,i,j) = _is_transition(a.data,basis(a),i,j)
function _is_transition(a::AbstractArray,b::NLevelBasis,i,j) 
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

function _is_transition(data,basis::CompositeBasis)
    N = length(basis.shape)
    ind = zeros(N)
    for k = 1:N
        if isa(basis.bases[k], NLevelBasis)
            nlevel_basis = basis.bases[k]
            # Check all possible transitions in the NLevelBasis
            for i in 1:nlevel_basis.N
                for j in 1:nlevel_basis.N
                    if i != j
                        ind .= 0
                        ind[k] = 1
                        transition_op = tensor([ (m!=k) ? identityoperator(basis.bases[m]) : transition(basis.bases[m],i,j) for (m,idx) in enumerate(ind)]...).data
                        
                        # Check if this is the same transition operator
                        if isequal(data, transition_op)
                            return k, i, j,1.0
                        end
                        
                        # Check if it's a scaled transition operator
                        non_zero_data = findall(!iszero, data)
                        non_zero_trans = findall(!iszero, transition_op)
                        
                        if non_zero_data == non_zero_trans && !isempty(non_zero_data)
                            scale_factor = data[first(non_zero_data)] / transition_op[first(non_zero_data)]
                            if all(idx -> isapprox(data[idx], transition_op[idx] * scale_factor), non_zero_data)
                                return k, i, j, scale_factor
                            end
                        end
                    end
                end
            end
        end
    end
    return 0, 0, 0,1.0
end

# Then modify each tensor function to include NLevelBasis checks



function tensor(a::T,b::Operator) where {T<:WaveguideOperatorT}
    if _is_identity(b)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,[a.indices...],(a.operators...,),a.factor)
    else
        # Check if b is a destroy operator
        k, scale_factor = _is_destroy(b.data,basis(b))
        if k > 0
            emission(a,basis(b),k;factor=scale_factor)
        else
            # Check if b is a create operator
            k, scale_factor = _is_create(b.data,basis(b))
            if k > 0
                absorption(a,basis(b),k;factor=scale_factor)
            else
                # Check for transition operators
                k, i, j, scale_factor = _is_transition(b.data,basis(b))
                if k > 0
                    btotal = basis(a) ⊗ basis(b)
                    return NLevelWaveguideOperator(btotal,btotal,a.factor*scale_factor,a.operators[1],[a.indices[1],k+length(basis(a).shape)],i,j)
                elseif !(typeof(basis(b)) <: QuantumOpticsBase.CompositeBasis)
                    btotal = basis(a) ⊗ basis(b)
                    LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
                else
                    error("I was trying to do LazyTensor(b.basis_l,b.basis_r,[1],(b,),1) the operator b in tensor(a::WaveguideOperator,b) seems to already be a tensor product (b.basis_l is a composite basis). Try reconstructing your tensor product so that WaveguideOperator is multiplied with other operators as the first step. Alternatively create the operator with LazyTensor(btotal,btotal,(1,2,3,...),(a,b,c,d...)) where btotal = b1 ⊗ b2 ⊗ b3 ⊗ ...")
                end
            end
        end
    end
end
function tensor(a::Operator,b::T) where {T<:WaveguideOperatorT}
    if _is_identity(a)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,[b.indices...].+1 ,(b.operators...,),b.factor)
    else
        # Check if a is a destroy operator
        k, scale_factor = _is_destroy(a.data,basis(a))
        if k > 0
            emission(basis(a),b,k;factor=scale_factor)
        else
            # Check if a is a create operator
            k, scale_factor = _is_create(a.data,basis(a))
            if k > 0
                absorption(basis(a),b,k;factor=scale_factor)
            else
                # Check for transition operators
                k, i, j, scale_factor = _is_transition(a.data,basis(a))
                if k > 0
                    btotal = basis(a) ⊗ basis(b)
                    return NLevelWaveguideOperator(btotal,btotal,b.factor*scale_factor,b.operators[1],[b.indices[1]+1,k+length(basis(b).shape)],i,j)
                elseif !(typeof(basis(a)) <: QuantumOpticsBase.CompositeBasis)
                    btotal = basis(a) ⊗ basis(b)
                    LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
                else
                    error("I was trying to do LazyTensor(b.basis_l,b.basis_r,[1],(b,),1) the operator b in tensor(a::WaveguideOperator,b) seems to already be a tensor product (b.basis_l is a composite basis). Try reconstructing your tensor product so that WaveguideOperator is multiplied with other operators as the first step. Alternatively create the operator with LazyTensor(btotal,btotal,(1,2,3,...),(a,b,c,d...)) where btotal = b1 ⊗ b2 ⊗ b3 ⊗ ...")
                end
            end
        end
    end
end

function tensor(a::T,b::Operator)  where {T<:WaveguideOperator}
    if _is_identity(b)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,(1,),(a,))
    else
        # Check if b is a destroy operator
        k, scale_factor = _is_destroy(b.data,basis(b))
        if k > 0
            emission(a,basis(b),k;factor=scale_factor)
        else
            # Check if b is a create operator
            k, scale_factor = _is_create(b.data,basis(b))
            if k > 0
                absorption(a,basis(b),k;factor=scale_factor)
            else
                # Check for transition operators
                k, i, j, scale_factor = _is_transition(b.data,basis(b))
                if k > 0
                    btotal = basis(a) ⊗ basis(b)
                    return NLevelWaveguideOperator(btotal,btotal,complex(scale_factor),a,[1,k+1],i,j)
                elseif !(typeof(basis(b)) <: QuantumOpticsBase.CompositeBasis)
                    btotal = basis(a) ⊗ basis(b)
                    LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
                else
                    error("I was trying to do LazyTensor(b.basis_l,b.basis_r,[1],(b,),1) the operator b in tensor(a::WaveguideOperator,b) seems to already be a tensor product (b.basis_l is a composite basis). Try reconstructing your tensor product so that WaveguideOperator is multiplied with other operators as the first step. Alternatively create the operator with LazyTensor(btotal,btotal,(1,2,3,...),(a,b,c,d...)) where btotal = b1 ⊗ b2 ⊗ b3 ⊗ ...")
                end
            end
        end
    end
end
function tensor(a::Operator,b::T) where {T<:WaveguideOperator}
    if _is_identity(a)
        btotal = basis(a) ⊗ basis(b)
        return LazyTensor(btotal,btotal,(length(basis(a).shape)+1,),(b,))
    else
        # Check if a is a destroy operator
        k, scale_factor = _is_destroy(a.data,basis(a))
        if k > 0
            emission(basis(a),b,k;factor=scale_factor)
        else
            # Check if a is a create operator
            k, scale_factor = _is_create(a.data,basis(a))
            if k > 0
                absorption(basis(a),b,k;factor=scale_factor)
            else
                # Check for transition operators
                k, i, j, scale_factor = _is_transition(a.data,basis(a))
                if k > 0
                    btotal = basis(a) ⊗ basis(b)
                    return NLevelWaveguideOperator(btotal,btotal,complex(scale_factor),b,[get_waveguide_location(basis(b))+length(basis(a).shape),k],i,j)
                elseif !(typeof(basis(a)) <: QuantumOpticsBase.CompositeBasis)
                    btotal = basis(a) ⊗ basis(b)
                    LazyTensor(btotal,btotal,(1,length(basis(a).shape)+1),(a,b))
                else
                    error("I was trying to do LazyTensor(b.basis_l,b.basis_r,[1],(b,),1) the operator b in tensor(a::WaveguideOperator,b) seems to already be a tensor product (b.basis_l is a composite basis). Try reconstructing your tensor product so that WaveguideOperator is multiplied with other operators as the first step. Alternatively create the operator with LazyTensor(btotal,btotal,(1,2,3,...),(a,b,c,d...)) where btotal = b1 ⊗ b2 ⊗ b3 ⊗ ...")
                end
            end
        end
    end
end


#TO DO: CLEAN UP
#Currently indexing is very complicated using the CavityIndexing structure and could possible be done smoother and faster.
"""
    mul!(result::Ket{B1}, a::CavityWaveguideEmission, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis}
    mul!(result::Ket{B1}, a::CavityWaveguideAbsorption, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis}
    
Fast in-place multiplication of operators/state vectors. Updates `result` as `result = alpha*a*b + beta*result`. `a` is a [`CavityWaveguideOperator`](@ref).
"""
function mul!(result::Ket{B1,A1}, a::CavityWaveguideEmission, b::Ket{B2,A2}, alpha, beta) where {B1<:Basis,B2<:Basis, A1<:AbstractArray, A2<:AbstractArray}
    dims = basis(result).shape
    i,j = a.loc[1],a.loc[2]
    beta1 = beta
    if iszero(beta)
        fill!(result.data,beta)
        beta1 = one(beta)
    end
    if length(dims) == 2
        loop_destroy_ax!(result.data,a,b.data,alpha,beta1,a.indexing)                     
    elseif length(dims) == 3
        loop_destroy_third_axis!(result.data,a,b.data,alpha,beta1,a.indexing,1)     
        a.indexing.idx_vec1[j] = dims[j]
        loop_rmul_axis!(result.data,a,b.data,alpha,beta1,a.indexing,1)
    else
        iterate_over_iter!(result.data,a,b.data,alpha,beta1,a.indexing,1,loop_destroy_third_axis!)
        a.indexing.idx_vec1[j] = dims[j]
        iterate_over_iter!(result.data,a,b.data,alpha,beta1,a.indexing,1,loop_rmul_axis!)
    end
    return result
end
function mul!(result::Ket{B1,A1}, a::CavityWaveguideAbsorption, b::Ket{B2,A2}, alpha, beta) where {B1<:Basis,B2<:Basis, A1<:AbstractArray, A2<:AbstractArray}
    dims = basis(result).shape
    i,j = a.loc[1],a.loc[2]
    beta1 = beta
    if iszero(beta)
        fill!(result.data,beta)
        beta1 = one(beta)
    end
    if length(dims) == 2
        loop_create_ax!(result.data,a,b.data,alpha,beta1,a.indexing)      
    elseif length(dims) == 3
        loop_create_third_axis!(result.data,a,b.data,alpha,beta1,a.indexing,1)
        a.indexing.idx_vec1[j] = 1
    else
        iterate_over_iter!(result.data,a,b.data,alpha,beta1,a.indexing,1,loop_create_third_axis!)
        a.indexing.idx_vec1[j] = 1
        iterate_over_iter!(result.data,a,b.data,alpha,beta1,a.indexing,1,loop_rmul_axis!)
    end
    return result
end

function loop_create_ax!(result,a,b,alpha,beta,indexing::WaveguideIndexing)
    i,j = a.loc[1],a.loc[2]
    for k in 2:indexing.k_idx
        indexing.idx_vec1[j] = k
        indexing.idx_vec2[j] = k-1
        li1 = linear_index(indexing.ndims,indexing.idx_vec1,indexing.strides)
        li2 = linear_index(indexing.ndims,indexing.idx_vec2,indexing.strides)  
        waveguide_mul!(view(result,li1:indexing.strides[a.loc[1]]:li1+(indexing.end_idx-1)*indexing.strides[a.loc[1]]),a.op,view(b,li2:indexing.strides[a.loc[1]]:li2+(indexing.end_idx-1)*indexing.strides[a.loc[1]]),sqrt(k-1)*a.factor*alpha,beta)
    end
end
function loop_destroy_ax!(result,a,b,alpha,beta,indexing::WaveguideIndexing)
    i,j = a.loc[1],a.loc[2]
    for k in 1:indexing.k_idx-1
        indexing.idx_vec1[j] = k
        indexing.idx_vec2[j] = k+1
        li1 = linear_index(indexing.ndims,indexing.idx_vec1,indexing.strides)
        li2 = linear_index(indexing.ndims,indexing.idx_vec2,indexing.strides)  
        waveguide_mul!(view(result,li1:indexing.strides[a.loc[1]]:li1+(indexing.end_idx-1)*indexing.strides[a.loc[1]]),a.op,view(b,li2:indexing.strides[a.loc[1]]:li2+(indexing.end_idx-1)*indexing.strides[a.loc[1]]),sqrt(k)*a.factor*alpha,beta)
    end
end
function loop_destroy_third_axis!(result,a,b,alpha,beta,indexing::WaveguideIndexing,idx)
    for j in 1:indexing.iter[idx]
        @inbounds indexing.idx_vec1[indexing.range_idx[idx]] = j
        @inbounds indexing.idx_vec2[indexing.range_idx[idx]] = j
        loop_destroy_ax!(result,a,b,alpha,beta,indexing)
    end
end
function loop_create_third_axis!(result,a,b,alpha,beta,indexing::WaveguideIndexing,idx)
    for j in 1:indexing.iter[idx]
        @inbounds indexing.idx_vec1[indexing.range_idx[idx]] = j
        @inbounds indexing.idx_vec2[indexing.range_idx[idx]] = j
        loop_create_ax!(result,a,b,alpha,beta,indexing)
    end
end
function loop_rmul_axis!(result,a,b,alpha,beta,indexing::WaveguideIndexing,idx)
    for j in 1:indexing.iter[idx]
        @inbounds indexing.idx_vec1[indexing.range_idx[idx]] = j
        li1 = linear_index(indexing.ndims,indexing.idx_vec1,indexing.strides)  
        rmul!(view(result,li1:a.indexing.strides[a.loc[1]]:li1+(a.indexing.end_idx-1)*a.indexing.strides[a.loc[1]]),beta)            
    end
end

function iterate_over_iter!(result,b,a,alpha,beta,indexing::WaveguideIndexing,idx,f::Function)
    if length(indexing.iter) == idx
        f(result,b,a,alpha,beta,indexing,idx)
    else
    for j in 1:indexing.iter[idx]
        @inbounds indexing.idx_vec1[indexing.range_idx[idx]] = j
        @inbounds indexing.idx_vec2[indexing.range_idx[idx]] = j
        iterate_over_iter!(result,b,a,alpha,beta,indexing,idx+1,f)
    end
    end
end

function iterate_rmul!(result,beta,iter::Tuple,idx_vec1,range_idx,idx,stride,end_idx,dims)
    if length(iter) == idx
        for i in Base.OneTo(iter[idx])
            idx_vec1[range_idx[idx]] = i
            li1 = linear_index(dims,idx_vec1)
            rmul!(view(result,li1:stride:li1+(end_idx-1)*stride),beta)            
        end
    else
        for j in Base.OneTo(iter[idx])
            idx_vec1[range_idx[idx]] = j
            iterate_rmul!(result,beta,iter,idx_vec1,range_idx,idx+1,stride,end_idx,dims)
        end
    end
end

function linear_index(dims::Vector, indices::Vector)
    n = length(dims)
    prev_prod = 1
    linear_idx = indices[1]
    for i = 2:n
        @inbounds prev_prod *= dims[i-1]
        @inbounds linear_idx += (indices[i] - 1) * prev_prod
    end
    return linear_idx
end

function linear_index(n::Int64, indices::Vector,strides::Vector)
    linear_idx = indices[1]
    for i = 2:n
        @inbounds linear_idx += (indices[i] - 1) * strides[i]
    end
    return linear_idx
end

function get_waveguide_operators(op::CavityWaveguideOperator)
    [op.op]
end
function set_waveguidetimeindex!(op::CavityWaveguideOperator,timeindex)
    op.op.timeindex = timeindex
end


#The following is experimental code to do indexing without reshaping.
"""
function mul!(result::Bra{B1}, a::Bra{B2}, b::CavityWaveguideOperator, alpha, beta) where {B1<:Basis,B2<:Basis}
    if a.loc[1] == 1 && a.loc[2] == 2
        return matmul_first!(result.data, b, a.data, alpha, beta)
    elseif a.loc[1] == length(b.basis.shape) && a.loc[2] == length(b.basis.shape)-1
        return matmul_last!(result.data, b, a.data, alpha, beta)
    end
    error("Waveguide operators in WaveguideQED operators has to be first or last")
end

   
#Calls CavityWaveguideoperator  on correct view of subsets of the state.
function matmul_first!(result, a::CavityWaveguideOperator, b, α::Number, β::Number)
    basis_shape = QuantumOpticsBase._comp_size(a.basis_l)
    d_first = basis_shape[1]*basis_shape[2]
    d_rest = length(b)÷d_first
    @uviews result b begin  # avoid allocations on reshape (NOT SURE IF WORKING)
        br = reshape(b, (d_first, d_rest))
        result_r = reshape(result, (d_first, d_rest))
        apply_first_op!(result_r,a,br,α,β)
    end
    result
end

#Same as matmul_first! But indexed in another way.
function matmul_last!(result, a::CavityWaveguideOperator, b, α::Number, β::Number)
    basis_shape = QuantumOpticsBase._comp_size(a.basis_l)
    d_last = basis_shape[end]*basis_shape[end-1]
    d_rest = length(b)÷d_last
    @uviews result b begin  # avoid allocations on reshape (NOT SURE IF WORKING)
        br = reshape(b, (d_rest, d_last))
        result_r = reshape(result, (d_rest, d_last))
        apply_last_op!(result_r,a,br,α,β)
    end
    result
end

#Same as matmul_first! But indexed in another way.
function matmul_last_mid!(result, a::CavityWaveguideOperator, b, α::Number, β::Number)
    basis_shape = size(b)
    d_last = basis_shape[end]*basis_shape[end-1]
    d_rest = length(b)÷d_last
    bp = b.parent
    rp = result.parent
    @uviews bp rp begin  # avoid allocations on reshape (NOT SURE IF WORKING)
        br = reshape(bp, (d_rest, d_last))
        result_r = reshape(rp, (d_rest, d_last))
        apply_last_op!(result_r,a,br,α,β)
    end
    result
end

#Calls CavityWaveguideoperator  on correct view of subsets of the state.
function matmul_first_mid!(result, a::CavityWaveguideOperator, b, α::Number, β::Number)
    basis_shape = size(b)
    d_first = basis_shape[1]*basis_shape[2]
    d_rest = length(b)÷d_first
    bp = b.parent
    rp = result.parent
    @uviews bp rp begin  # avoid allocations on reshape (NOT SURE IF WORKING)
        br = reshape(bp, (d_first, d_rest))
        result_r = reshape(result, (d_first, d_rest))
        apply_first_op!(result_r,a,br,α,β)
    end
    result
end


function apply_last_op!(result,a::CavityWaveguideOperator,br,α,β)
    @simd for i in axes(br,1)
        waveguide_mul_last!(view(result,i,:), a, view(br,i,:), α, β)
    end
end

function apply_first_op!(result,a::CavityWaveguideOperator,br,α,β)
    @simd for i in axes(br,2)
        waveguide_mul_first!(view(result,:,i), a, view(br,:,i), α, β)
    end
end

function waveguide_mul_last!(result, a::CavityWaveguideAbsorption, b, alpha::Number, beta::Number)
    b_data = reshape(b, :, length(basis(a).bases[a.loc[1]]))
    result_data = reshape(result,  :, length(basis(a).bases[a.loc[1]]))
    @simd for i in 2:size(b_data,1)
        waveguide_mul!(view(result_data,i,:),a.op,view(b_data,i-1,:),sqrt(i-1)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,1,:),beta)
end

function waveguide_mul_first!(result, a::CavityWaveguideAbsorption, b, alpha::Number, beta::Number)
    b_data = reshape(b, length(basis(a).bases[a.loc[a.loc[1]]]),:)
    result_data = reshape(result,length(basis(a).bases[a.loc[1]]),:)
    @simd for i in 2:size(b_data,2)
        waveguide_mul!(view(result_data,:,i),a.op,view(b_data,:,i-1),sqrt(i-1)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,:,1),beta)
end

function waveguide_mul_last!(result, a::CavityWaveguideEmission, b, alpha::Number, beta::Number)
    b_data = reshape(b, :, length(basis(a).bases[a.loc[1]]))
    result_data = reshape(result,  :, length(basis(a).bases[a.loc[1]]))
    for i in 1:size(b_data,1)-1
        waveguide_mul!(view(result_data,i,:),a.op,view(b_data,i+1,:),sqrt(i)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,size(b_data,1),:),beta)
end

function waveguide_mul_first!(result, a::CavityWaveguideEmission, b, alpha::Number, beta::Number)
    b_data = reshape(b,length(basis(a).bases[a.loc[1]]),:)
    result_data = reshape(result,length(basis(a).bases[a.loc[1]]),:)
    for i in 1:size(b_data,2)-1
        waveguide_mul!(view(result_data,:,i),a.op,view(b_data,:,i+1),sqrt(i)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,:,size(b_data,2)),beta)
end


function matmul_first!(result, a::CavityWaveguideOperator, b, α::Number, β::Number)
    basis_shape = QuantumOpticsBase._comp_size(a.basis_l)
    d_first = basis_shape[1]*basis_shape[2]
    d_rest = length(b)÷d_first
    mapping = reshape(1:1:length(result), (d_first, d_rest))
    apply_first_op!(result,a,b,α,β,mapping)
    result
end

"""

"""
function matmul_last!(result, a::CavityWaveguideOperator, b, α::Number, β::Number)
    basis_shape = QuantumOpticsBase._comp_size(a.basis_l)
    d_first = basis_shape[end]*basis_shape[end-1]
    d_rest = length(b)÷d_first
    mapping = reshape(1:1:length(result), (d_rest, d_first))
    apply_last_op!(result,a,b,α,β,mapping)
    result
end


function apply_last_op!(result,a::CavityWaveguideOperator,br,α,β,mapping)
    @simd for i in axes(mapping,1)
        waveguide_mul_last!(result, a, br, α, β,view(mapping,i,:))
    end
end


function apply_first_op!(result,a::CavityWaveguideOperator,br,α,β,mapping)
    for i in axes(mapping,2)
        waveguide_mul_first!(result, a, br, α, β,view(mapping,:,i))
    end
end


function waveguide_mul_first!(result, a::CavityWaveguideAbsorption, b, alpha::Number, beta::Number,mapping)
    mapping = reshape(mapping,length(basis(a).bases[1]),:)
    for i in 2:size(mapping,2)
        waveguide_mul!(view(result,mapping[1,i]:mapping[2,i]-mapping[1,i]:mapping[end,i]),a.op,view(b,mapping[1,i-1]:mapping[2,i-1]-mapping[1,i-1]:mapping[end,i-1]),sqrt(i-1)*alpha*a.factor,beta)
    end
    rmul!(view(result,mapping[1,1]:mapping[2,1]-mapping[1,1]:mapping[end,1]),beta)
end


function waveguide_mul_last!(result, a::CavityWaveguideEmission, b, alpha::Number, beta::Number,mapping)
    map = reshape(mapping,:,length(basis(a).bases[end]))
    @simd for i in 1:size(map,1)-1
        waveguide_mul!(view(result,map[i,1]:map[i,2]-map[i,1]:map[i,end]),a.op,view(b,map[i+1,1]:map[i+1,2]-map[i+1,1]:map[i+1,end]),sqrt(i-1)*alpha*a.factor,beta)
    end
    rmul!(view(result,map[1,1]:map[1,2]-map[1,1]:map[1,end]),beta)
end


function waveguide_mul_first!(result, a::CavityWaveguideEmission, b, alpha::Number, beta::Number,mapping)
    map = reshape(mapping,length(basis(a).bases[1]),:)
    @simd for i in 1:size(map,2)-1
        waveguide_mul!(view(result,map[1,i]:map[2,i]-map[1,i]:map[end,i]),a.op,view(b,map[1,i+1]:map[2,i+1]-map[1,i+1]:map[end,i+1]),sqrt(i-1)*alpha*a.factor,beta)
    end
    rmul!(view(result,map[1,1]:map[1,2]-map[1,1]:map[1,end]),beta)
end
"""