"""
Abstract class for WaveguideOperators. Used to dispatch special mul! function.
"""
abstract type WaveguideOperator{B1,B2} <: AbstractOperator{B1,B2} end


"""
    WaveguideDestroy{B1,B2,Np,idx} <: WaveguideOperator{B1,B2}

Operator structure for dispatching annihilation operation on Waveguide state.
Np is used to dispatch one or two photon routine and idx denotes the index of the waveguide the operator is acting on. 
"""
mutable struct WaveguideDestroy{B1,B2,Np,idx} <: WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end

"""
    WaveguideCreate{B1,B2,N,idx} <: WaveguideOperator{B1,B2}

Operator structure for dispatching creation operation on Waveguide state.
Np is used to dispatch one or two photon routine and idx denotes the index of the waveguide the operator is acting on. 
"""
mutable struct WaveguideCreate{B1,B2,N,idx} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end


function Base.:eltype(x::WaveguideOperator) typeof(x.factor) end

#Base.:*(x::WaveguideOperator{B1,B2},y::WaveguideOperator{B1,B2}) where {B1,B2} = LazyProduct((x,y),x.factor*y.factor)


#Methods for copying waveguide operators
function Base.:copy(x::WaveguideDestroy{B,B,Np,idx}) where {B,Np,idx}
    WaveguideDestroy{B,B,Np,idx}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end
function Base.:copy(x::WaveguideCreate{B,B,Np,idx}) where {B,Np,idx}
    WaveguideCreate{B,B,Np,idx}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end

#Arithmetic operations for multiplying, which updates factor in the operator.
function Base.:*(a::Number,b::WaveguideOperator)
    out = copy(b)
    out.factor=out.factor*a
    return out
end
Base.:*(b::WaveguideOperator,a::Number)=*(a,b)
function Base.:/(a::WaveguideOperator,b::Number)
    out = copy(a)
    out.factor=out.factor/b
    return out
end
Base.:-(a::WaveguideOperator) = *(-1,a)

"""
    destroy(basis::WaveguideBasis{Np,1}) where {Np}
    destroy(basis::WaveguideBasis{Np,Nw},i::Int) where {Np,Nw}

Annihilation operator ``w`` for [`WaveguideBasis`](@ref) ``w_i(t_k) | 1_j \\emptyset \\rangle_i = \\delta_{k,j} | \\emptyset \\rangle`` with cutoff (maximum number of photons Np) where `i` is the index of the waveguide.
``t_k`` is determined by the timeindex property of the operator which can be changed by [`set_waveguidetimeindex!(a::WaveguideOperator,k::Int)`](@ref) 

# Arguments
- basis of type WaveguideBasis, defines cutoff photon number `Np` and number of waveguides `Nw`
- `i` determines which waveguide the operator acts on and should be `i ≤ Nw`. If `Nw=1` then `i=1` is assumed (there is only one waveguide).

# Returns
[`WaveguideDestroy`](@ref)
"""
function destroy(basis::WaveguideBasis{Np,1}) where {Np}
    B = typeof(basis)
    return WaveguideDestroy{B,B,Np,1}(basis,basis,1,1)
end
function destroy(basis::WaveguideBasis{Np,Nw},i::Int) where {Np,Nw}
    @assert i <= Nw
    B = typeof(basis)
    return WaveguideDestroy{B,B,Np,i}(basis,basis,1,1)
end

"""
    create(basis::WaveguideBasis{Np,1}) where {Np}
    create(basis::WaveguideBasis{Np,Nw},i::Int) where {Np,Nw}

Creation operator ``w^\\dagger`` for [`WaveguideBasis`](@ref) ``w_i(t_k)^\\dagger | \\emptyset \\rangle = | 1_k \\emptyset \\rangle_i `` with cutoff (maximum number of photons Np) where `i` is the index of the waveguide.
``t_k`` is determined by the timeindex property of the operator which can be changed by [`set_waveguidetimeindex!(a::WaveguideOperator,k::Int)`](@ref) 


# Arguments
- basis of type WaveguideBasis, defines cutoff photon number `Np` and number of waveguides `Nw`
- `i` determines which waveguide the operator acts on and should be `i ≤ Nw`. If `Nw=1` then `i=1` is assumed (there is only one waveguide).

# Returns
[`WaveguideCreate`](@ref)
"""
function create(basis::WaveguideBasis{Np,1}) where {Np}
    B = typeof(basis)
    return WaveguideCreate{B,B,Np,1}(basis,basis,1,1)
end
function create(basis::WaveguideBasis{Np,Nw},i::Int) where {Np,Nw}
    @assert i <= Nw
    B = typeof(basis)
    return WaveguideCreate{B,B,Np,i}(basis,basis,1,1)
end


"""
    dagger(op::WaveguideCreate)
    dagger(op::WaveguideCreate)

Dagger opration on Waveguide operator. 
"""
function dagger(op::WaveguideCreate{B,B,Np,idx}) where {B,Np,idx}
    @assert op.basis_l == op.basis_r
    WaveguideDestroy{B,B,Np,idx}(op.basis_l,op.basis_r,op.factor,op.timeindex)
end
function dagger(op::WaveguideDestroy{B,B,Np,idx}) where {B,Np,idx}
    @assert op.basis_l == op.basis_r
    WaveguideCreate{B,B,Np,idx}(op.basis_l,op.basis_r,op.factor,op.timeindex)
end


"""
    tensor(op1::AbstractOperator,op2::WaveguideOperator)
    tensor(op1::WaveguideOperator,op2::AbstractOperator)

Methods for tensorproducts between QuantumOptics.jl operator and [`WaveguideOperator`](@ref). This is done by forming a LazyTensor.
"""
#TODO: Update method to allow for three or more hilbert spaces. Should be fixed
function tensor(a::DataOperator,b::WaveguideOperator) 
    if isequal(a,identityoperator(basis(a)))
        btotal = basis(a) ⊗ basis(b)
        LazyTensor(btotal,btotal,[length(basis(a).shape)+1],(b,),1)
    else
        LazyTensor(a.basis_l,a.basis_r,[1],(a,),1) ⊗ LazyTensor(b.basis_l,b.basis_r,[1],(b,),1)
    end
end
function tensor(a::WaveguideOperator,b::DataOperator) 
    if isequal(b,identityoperator(basis(b)))
        btotal = basis(a) ⊗ basis(b)
        LazyTensor(btotal,btotal,[1],(a,),1)
    else
        LazyTensor(a.basis_l,a.basis_r,[1],(a,),1) ⊗ LazyTensor(b.basis_l,b.basis_r,[1],(b,),1)
    end
end
 
"""
    identityoperator(a::WaveguideOperator)

Return identityoperator(a.basis_l).
"""
function QuantumOpticsBase.identityoperator(a::WaveguideOperator)
    identityoperator(a.basis_l)
end
function QuantumOpticsBase.identityoperator(::Type{T}, b1::Basis, b2::Basis) where {T<:WaveguideOperator}
    @assert b1==b2
    identityoperator(b1)
end
function identityoperator(::Type{T},::Type{ComplexF64}, b1::Basis, b2::Basis) where {T<:WaveguideOperator}
    @assert b1==b2
    identityoperator(b1)
end

"""
    mul!(result::Ket{B1}, a::LazyTensor{B1,B2,F,I,T}, b::Ket{B2}, alpha, beta)
    
In-place multiplication of operators/state vectors. Updates `result` as `result = alpha*a*b + beta*result`. `a` is a LazyTensor that contains a [`WaveguideOperator`](@ref) 
"""
function mul!(result::Ket{B1}, a::LazyTensor{B1,B2,F,I,T}, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis, F,I,T<:Tuple{Vararg{AbstractOperator}}}
    b_data = Base.ReshapedArray(b.data, QuantumOpticsBase._comp_size(basis(b)), ())
    result_data = Base.ReshapedArray(result.data, QuantumOpticsBase._comp_size(basis(result)), ())

    tp_ops = QuantumOpticsBase._tpops_tuple(a)
    iso_ops = QuantumOpticsBase._explicit_isometries(eltype(a), a.indices, a.basis_l, a.basis_r)

    QuantumOpticsBase._tp_sum_matmul!(result_data, tp_ops, iso_ops, b_data, alpha * a.factor, beta)
    result
end

function QuantumOpticsBase._tpops_tuple(a::LazyTensor{B1,B2,F,I,T}; shift=0, op_transform=identity)  where {B1<:Basis,B2<:Basis, F,I,T<:Tuple{Vararg{AbstractOperator}}}
    length(a.operators) == 0 == length(a.indices) && return ()
    op_pairs = tuple((((isa(op,DataOperator) ? op_transform(op.data) : op_transform(op)), i + shift) for (op, i) in zip(a.operators, a.indices))...)

    # Filter out identities:
    # This induces a non-trivial cost only if _is_square_eye is not inferrable.
    # This happens if we have Eyes that are not SquareEyes.
    # This can happen if the user constructs LazyTensor operators including
    # explicit identityoperator(b,b).
    filtered = filter(p->!QuantumOpticsBase._is_square_eye(p[1]), op_pairs)
    return filtered
end
QuantumOpticsBase._is_square_eye(a::AbstractOperator) = false


#Called from _tp_sum_matmul!
#Makes sure operator works on correct part of tensor.
function QuantumOpticsBase.:_tp_matmul!(result, a::WaveguideOperator, loc::Integer, b, α::Number, β::Number)
    if loc == 1
        return QuantumOpticsBase._tp_matmul_first!(result, a, b, α, β)
    elseif loc == ndims(b)
        return QuantumOpticsBase._tp_matmul_last!(result, a, b, α, β)
    end
    QuantumOpticsBase._tp_matmul_mid!(result, a, loc, b, α, β)
end
function QuantumOpticsBase._tp_matmul_mid!(result, a::WaveguideOperator, loc::Integer, b, α::Number, β::Number)
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

#Called from _tp_matmul! and _tp_matmul_mid!
#Calls waveguide_mul! on correct view of whole subsets of the state.
function QuantumOpticsBase.:_tp_matmul_first!(result::Base.ReshapedArray, a::WaveguideOperator, b::Base.ReshapedArray, α::Number, β::Number)
    d_first = size(b, 1)
    d_rest = length(b)÷d_first
    bp = b.parent
    rp = result.parent
    @uviews bp rp begin  # avoid allocations on reshape
        br = reshape(bp, (d_first, d_rest))
        result_r = reshape(rp, (size(a, 1), d_rest))
        apply_first_op!(result_r,a,br,α,β)
    end
    result
end

#Same as _tp_matmul_first! But indexed in another way.
function QuantumOpticsBase.:_tp_matmul_last!(result::Base.ReshapedArray, a::WaveguideOperator, b::Base.ReshapedArray, α::Number, β::Number)
    d_last = size(b, ndims(b))
    d_rest = length(b)÷d_last
    bp = b.parent
    rp = result.parent
    @uviews bp rp begin  # avoid allocations on reshape
        br = reshape(bp, (d_rest, d_last))
        result_r = reshape(rp, (d_rest, size(a, 1)))
        apply_last_op!(result_r,a,br,α,β)
    end
    result
end

function apply_last_op!(result,a::WaveguideOperator,br,α,β)
    for i in axes(br,1)
        waveguide_mul!(view(result,i,:), a, view(br,i,:), α, β)
    end
end
function apply_first_op!(result,a::WaveguideOperator,br,α,β)
    for i in axes(br,2)
        waveguide_mul!(view(result,:,i), a, view(br,:,i), α, β)
    end
end

#Get tempory vector
function QuantumOpticsBase._tp_sum_get_tmp(op::WaveguideOperator, loc::Integer, arr::AbstractArray{S,N}, sym) where {S,N}
    shp = ntuple(i -> i == loc ? size(op,1) : size(arr,i), N)
    QuantumOpticsBase._tp_matmul_get_tmp(S, shp, sym,arr)
end


function mul!(result::Ket{B},a::WaveguideOperator{B,B},input::Ket{B},alpha,beta) where {B}
    waveguide_mul!(result.data,a,input.data,alpha,beta)
end

#Destroy 1 waveguide photon
function waveguide_mul!(result,a::WaveguideDestroy{B,B,1,idx},b,alpha,beta) where {B,idx}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    @inbounds result[1] += alpha*a.factor*b[a.timeindex+(idx-1)*a.basis_l.nsteps+1]
    return
end
#Destroy 2 waveguide photon
function waveguide_mul!(result,a::WaveguideDestroy{B,B,2,1},b,alpha,beta) where {B<:SingleWaveguideBasis}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += alpha*a.factor*b[timeindex+1]
    #twophotonview = TwoPhotonTimestepView(b,timeindex,nsteps,nsteps+1)
    #axpy!(alpha*a.factor,twophotonview,view(result,2:1:nsteps+1))
    twophoton_destroy!(view(result,2:1:nsteps+1),b,alpha*a.factor,timeindex,nsteps,nsteps+1)
    return
end
#Destroy 2 waveguide photon
function waveguide_mul!(result,a::WaveguideDestroy{B,B,2,idx},b,alpha,beta) where {B<:MultipleWaveguideBasis,idx}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw  = get_number_of_waveguides(a.basis_l)
    @inbounds result[1] += alpha*a.factor*b[timeindex+(idx-1)*a.basis_l.nsteps+1]
    #twophotonview = TwoPhotonTimestepView(b,timeindex,nsteps,Nw*nsteps+1+(idx-1)*(nsteps*(nsteps+1))÷2)
    #axpy!(alpha*a.factor,twophotonview,view(result,2+(idx-1)*nsteps:1:idx*nsteps+1))
    twophoton_destroy!(view(result,2+(idx-1)*nsteps:1:idx*nsteps+1),b,alpha*a.factor,timeindex,nsteps,Nw*nsteps+1+(idx-1)*(nsteps*(nsteps+1))÷2)
    for k in filter(x -> x != idx, 1:Nw)
        i,j = min(k,idx),max(k,idx)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        #twophotonview = TwoWaveguideTimestepView(b,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx)
        #axpy!(alpha*a.factor,twophotonview,view(result,2+(k-1)*nsteps:1:k*nsteps+1))
        twowaveguide_destroy!(view(result,2+(k-1)*nsteps:1:k*nsteps+1),b,alpha*a.factor,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx)
    end
    return
end


#Create 1 waveguide photon 
function waveguide_mul!(result,a::WaveguideCreate{B,B,1,idx},b,alpha,beta) where {B,idx}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    @inbounds result[a.timeindex+(idx-1)*a.basis_l.nsteps+1] += alpha*a.factor*b[1]
    return
end
#Create 2 waveguide photon
function waveguide_mul!(result,a::WaveguideCreate{B,B,2,1},b,alpha,beta) where {B<:SingleWaveguideBasis}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[timeindex+1] += alpha*a.factor*b[1]
    #twophotonview = TwoPhotonTimestepView(result,timeindex,nsteps,nsteps+1)
    #axpy!(alpha*a.factor,view(b,2:1:nsteps+1),twophotonview)
    twophoton_create!(result,view(b,2:1:nsteps+1),alpha*a.factor,timeindex,nsteps,nsteps+1)
    return
end
#Create 2 waveguide photon
function waveguide_mul!(result,a::WaveguideCreate{B,B,2,idx},b,alpha,beta) where {B<:MultipleWaveguideBasis,idx}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw  = get_number_of_waveguides(a.basis_l)
    @inbounds result[timeindex+(idx-1)*a.basis_l.nsteps+1] += alpha*a.factor*b[1]
    #twophotonview = TwoPhotonTimestepView(result,timeindex,nsteps,Nw*nsteps+1+(idx-1)*(nsteps*(nsteps+1))÷2)
    #axpy!(alpha*a.factor,view(b,2+(idx-1)*nsteps:1:idx*nsteps+1),twophotonview)
    twophoton_create!(result,view(b,2+(idx-1)*nsteps:1:idx*nsteps+1),alpha*a.factor,timeindex,nsteps,Nw*nsteps+1+(idx-1)*(nsteps*(nsteps+1))÷2)
    @simd for k in filter(x -> x != idx, 1:Nw)
        i,j = min(k,idx),max(k,idx)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        #twophotonview = TwoWaveguideTimestepView(result,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx)
        #axpy!(alpha*a.factor,view(b,2+(k-1)*nsteps:1:k*nsteps+1),twophotonview)
        twowaveguide_create!(result,view(b,2+(k-1)*nsteps:1:k*nsteps+1),alpha*a.factor,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx)
    end
    return
end

function twophoton_create!(result,b,alpha,timeindex,nsteps,offset)
    @simd for i in 1:timeindex-1
        @inbounds result[offset + twophoton_index(i,nsteps,timeindex)] += alpha*b[i] 
    end
    @simd for i in timeindex+1:lastindex(b)
        @inbounds result[offset + twophoton_index(timeindex,nsteps,i)] += alpha*b[i] 
    end
    @inbounds result[offset + twophoton_index(timeindex,nsteps,timeindex)] += sqrt(2)*alpha*b[timeindex]
end
function twophoton_destroy!(result,b,alpha,timeindex,nsteps,offset)
    @simd for i in 1:timeindex-1
        @inbounds result[i]  += alpha*b[offset + twophoton_index(i,nsteps,timeindex)]
    end
    @simd for i in timeindex+1:lastindex(result)
        @inbounds result[i]  += alpha*b[offset + twophoton_index(timeindex,nsteps,i)]
    end
    @inbounds result[timeindex] += sqrt(2)*alpha*b[offset + twophoton_index(timeindex,nsteps,timeindex)]
end

function twowaveguide_create!(result,b,alpha,timeindex,nsteps,offset,order)
    if order
        @simd for i in eachindex(b)
            @inbounds result[offset + (i-1)*nsteps + timeindex] += alpha*b[i] 
        end
    else
        @simd for i in eachindex(b)
            @inbounds result[offset + (timeindex-1)*nsteps + i] += alpha*b[i] 
        end
    end
end
function twowaveguide_destroy!(result,b,alpha,timeindex,nsteps,offset,order)
    if order
        @simd for i in eachindex(result)
            @inbounds result[i]  += alpha*b[offset + (i-1)*nsteps + timeindex]
        end
    else
        @simd for i in eachindex(result)
            @inbounds result[i]  += alpha*b[offset + (timeindex-1)*nsteps + i]
        end
    end
end





@inline function add_zerophoton_onephoton!(a,b,alpha,timeindex::Int)
    @inbounds a[1] += alpha*b[timeindex+1]
end

@inline function add_onephoton_zerophoton!(a,b,alpha,timeindex::Int)
    @inbounds a[1+timeindex] += alpha*b[1]
end

function add_twophoton_onephoton!(a,b,alpha)
    @simd for j in eachindex(a)
        @inbounds a[j] += alpha*b[j+1]
    end
end

function add_onephoton_twophoton!(a,b,alpha,nsteps::Int)
    @simd for j in 1:nsteps
        @inbounds a[j+1] += alpha*b[j]
    end
end


"""
    get_waveguide_operators(basis::LazySum)
    get_waveguide_operators(basis::LazyProduct)
    get_waveguide_operators(basis::LazyTensor)
    get_waveguide_operators(basis::Tuple)
    get_waveguide_operators(basis::Array)
    get_waveguide_operators(basis::WaveguideOperator)
    
Returns all [`WaveguideOperator`](@ref) in LazyOperator or from a list of operators. If no [`WaveguideOperator`](@ref) is found, and empty array is returned.

"""
function get_waveguide_operators(op::LazySum)
    out = [[if !isnothing(get_waveguide_operators(O)) get_waveguide_operators(O) end for O in op.operators]...]
    out = collect(Iterators.flatten(out))
    out[findall(x->typeof(x)<:WaveguideOperator,out)]
end
function get_waveguide_operators(op::LazyProduct)
    out = [[if !isnothing(get_waveguide_operators(O)) get_waveguide_operators(O) end for O in op.operators]...]
    out = collect(Iterators.flatten(out))
    out[findall(x->typeof(x)<:WaveguideOperator,out)]
end
function get_waveguide_operators(op::LazyTensor)
    out = [[if !isnothing(get_waveguide_operators(O)) get_waveguide_operators(O) end for O in op.operators]...]
    out = collect(Iterators.flatten(out))
    out[findall(x->typeof(x)<:WaveguideOperator,out)]
end
function get_waveguide_operators(op::Tuple)
    out = [[if !isnothing(get_waveguide_operators(O)) get_waveguide_operators(O) end for O in op]...]
    out = collect(Iterators.flatten(out))
    out[findall(x->typeof(x)<:WaveguideOperator,out)]
end
function get_waveguide_operators(op::Array)
    out = [[if !isnothing(get_waveguide_operators(O)) get_waveguide_operators(O) end for O in op]...]
    out = collect(Iterators.flatten(out))
    out[findall(x->typeof(x)<:WaveguideOperator,out)]
end
function get_waveguide_operators(op::WaveguideOperator)
    [op]
end
function get_waveguide_operators(op)
    []
end

"""
    get_waveguidetimeindex(op)

Return timeindex of operator or list of operators containing [`WaveguideOperator`](@ref) and assert that all timeindeces are the same. 
"""
function get_waveguidetimeindex(op)
    ops = get_waveguide_operators(op)
    timeindex = ops[1].timeindex
    for a in ops[2:end]
        @assert a.timeindex == timeindex
    end
    return timeindex
end
function get_waveguidetimeindex(op::WaveguideOperator)
    op.timeindex
end


"""
    set_waveguidetimeindex!(op,index)

Set timeindex of all [`WaveguideOperator`](@ref) in operator or list of operators to index
"""
function set_waveguidetimeindex!(op::Vector{T},index) where T<:WaveguideOperator
    for O in op
        set_waveguidetimeindex!(O,index)
    end
end
function set_waveguidetimeindex!(op::WaveguideOperator,index)
    op.timeindex = index
end
function set_waveguidetimeindex!(op::LazyProduct,index)
    for x in op.operators
        set_waveguidetimeindex!(x,index)
    end
end
function set_waveguidetimeindex!(op::LazySum,index)
    for x in op.operators
        set_waveguidetimeindex!(x,index)
    end
end
function set_waveguidetimeindex!(op::LazyTensor,index)
    for x in op.operators
        set_waveguidetimeindex!(x,index)
    end
end
function set_waveguidetimeindex!(op::Tuple,index)
    for x in op
        set_waveguidetimeindex!(x,index)
    end
end
function set_waveguidetimeindex!(op::Array,index)
    for x in op
        set_waveguidetimeindex!(x,index)
    end
end
function set_waveguidetimeindex!(args...)
    for i in 1:length(args)-1
        set_waveguidetimeindex!(args[i],args[end])
    end
end
function set_waveguidetimeindex!(op,index)

end
