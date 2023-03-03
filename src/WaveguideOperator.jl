"""
Abstract class for WaveguideOperators. Used to dispatch special mul! function.
"""
abstract type WaveguideOperator{B1,B2} <: AbstractOperator{B1,B2} end


"""
    WaveguideDestroy{N} <: WaveguideOperator

Operator structure for dispatching annihilation operation on Waveguide state.
N is used to dispatch one or two photon routine. 
"""
mutable struct WaveguideDestroy{B1,B2,N} <: WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end

"""
    WaveguideCreate{N} <: WaveguideOperator

Operator structure for dispatching creation operation on Waveguide state.
N is used to dispatch one or two photon routine. 
"""
mutable struct WaveguideCreate{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end


function Base.:eltype(x::WaveguideOperator) typeof(x.factor) end

#Methods for copying waveguide operators
function Base.:copy(x::WaveguideDestroy{B,B,1}) where {B}
    WaveguideDestroy{B,B,1}(x.basis_l,x.basis_r,x.factor,1)
end
function Base.:copy(x::WaveguideDestroy{B,B,2}) where {B}
    WaveguideDestroy{B,B,2}(x.basis_l,x.basis_r,x.factor,1)
end
function Base.:copy(x::WaveguideCreate{B,B,1}) where {B}
    WaveguideCreate{B,B,1}(x.basis_l,x.basis_r,x.factor,1)
end
function Base.:copy(x::WaveguideCreate{B,B,2}) where {B}
    WaveguideCreate{B,B,2}(x.basis_l,x.basis_r,x.factor,1)
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
    destroy(basis::WaveguideBasis{1})
    destroy(basis::WaveguideBasis{2})

Annihilation operator for [`WaveguideBasis`](@ref) for either one or two photons. 

"""
function destroy(basis::WaveguideBasis{1})
    B = typeof(basis)
    return WaveguideDestroy{B,B,1}(basis,basis,1,1)
end
function destroy(basis::WaveguideBasis{2})
    B = typeof(basis)
    return WaveguideDestroy{B,B,2}(basis,basis,1,1)
end

"""
    create(basis::WaveguideBasis{1})
    create(basis::WaveguideBasis{2})

Creation operator for [`WaveguideBasis`](@ref) for either one or two photons. 

"""
function create(basis::WaveguideBasis{1})
    B = typeof(basis)
    return WaveguideCreate{B,B,1}(basis,basis,1,1)
end
function create(basis::WaveguideBasis{2})
    B = typeof(basis)
    return WaveguideCreate{B,B,2}(basis,basis,1,1)
end

"""
    empty(basis:WaveguideBasis{1})
    empty(basis:WaveguideBasis{2})

Empty operator for [`WaveguideBasis`](@ref) for either one or two photons.
"""
function projector(basis::WaveguideBasis{1})
    return WaveguideProject{1}(basis,basis,1,1)
end
function projector(basis::WaveguideBasis{2})
    return WaveguideProject{2}(basis,basis,1,1)
end


"""
    dagger(basis::WaveguideBasis{1})
    dagger(basis::WaveguideBasis{2})

Dagger opration on Waveguide operator. 

"""
function dagger(op::WaveguideCreate)
    @assert op.basis_l == op.basis_r
    out = destroy(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out 
end
function dagger(op::WaveguideDestroy)
    @assert op.basis_l == op.basis_r
    out = create(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out
end


"""
    tensor(op1::AbstractOperator,op2::WaveguideOperator)
    tensor(op1::WaveguideOperator,op2::AbstractOperator)

Methods for tensorproducts between QuantumOptics.jl operator and [`WaveguideOperator`](@ref). This is done by forming a LazyTensor.
"""
#TODO: Update method to allow for three or more hilbert spaces.
function tensor(a::DataOperator,b::WaveguideOperator) 
    LazyTensor(a.basis_l,a.basis_r,[1],(a,),1) ⊗ LazyTensor(b.basis_l,b.basis_r,[1],(b,),1)
end
function tensor(a::WaveguideOperator,b::DataOperator) 
    LazyTensor(a.basis_l,a.basis_r,[1],(a,),1) ⊗ LazyTensor(b.basis_l,b.basis_r,[1],(b,),1)
end
 
"""
    identityoperator(a::WaveguideOperator)

Return identityoperator(a.basis_l).
QUESTION: (does basis_l or basis_r matter?)
"""
function QuantumOptics.identityoperator(a::WaveguideOperator)
    identityoperator(a.basis_l)
end
function QuantumOptics.identityoperator(::Type{T}, b1::Basis, b2::Basis) where {T<:WaveguideOperator}
    @assert b1==b2
    identityoperator(b1)
end

"""
    mul!(result::Ket{B1}, a::LazyTensor{B1,B2,F,I,T}, b::Ket{B2}, alpha, beta)
    mul!(result::Bra{B1}, a::Bra{B2}, b::LazyTensor{B1,B2,F,I,T}, alpha, beta)

In-place multiplication of operators/state vectors. Updates `result` as `result = alpha*a*b + beta*result`. `a` is a LazyTensor that contains a [`WaveguideOperator`](@ref) 

"""
function mul!(result::Ket{B1}, a::LazyTensor{B1,B2,F,I,T}, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis, F,I,T<:Tuple{Vararg{AbstractOperator}}}
    b_data = Base.ReshapedArray(b.data, QuantumOpticsBase._comp_size(basis(b)), ())
    result_data = Base.ReshapedArray(result.data, QuantumOpticsBase._comp_size(basis(result)), ())

    tp_ops = (tuple(( (isa(op,DataOperator) ? op.data : op) for op in a.operators)...), a.indices)
    iso_ops = QuantumOpticsBase._explicit_isometries(a.indices, a.basis_l, a.basis_r)

    QuantumOpticsBase._tp_sum_matmul!(result_data, tp_ops, iso_ops, b_data, alpha * a.factor, beta)
    result
end
function mul!(result::Bra{B1}, a::Bra{B2}, b::LazyTensor{B1,B2,F,I,T}, alpha, beta) where {B1<:Basis,B2<:Basis, F,I,T<:Tuple{Vararg{AbstractOperator}}}
    a_data = Base.ReshapedArray(a.data, QuantumOpticsBase._comp_size(basis(a)), ())
    result_data = Base.ReshapedArray(result.data, QuantumOpticsBase._comp_size(basis(result)), ())

    tp_ops = (tuple(( (isa(op,DataOperator) ? op.data : op) for op in b.operators)...), b.indices)
    iso_ops = QuantumOpticsBase._explicit_isometries(b.indices, b.basis_l, b.basis_r)

    QuantumOpticsBase._tp_sum_matmul!(result_data, tp_ops, iso_ops, a_data, alpha * b.factor, beta)
    result
end

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
    QuantumOpticsBase._tp_matmul_get_tmp(S, shp, sym)
end


function mul!(result::Ket{B},a::WaveguideOperator{B,B},input::Ket{B},alpha,beta) where {B}
    waveguide_mul!(result.data,a,input.data,alpha,beta)
end

#Destroy 1 waveguide photon
function waveguide_mul!(result,a::WaveguideDestroy{B,B,1},b,alpha,beta) where {B}
    rmul!(result,beta)
    add_zerophoton_onephoton!(result,b,alpha*a.factor,a.timeindex)
    return
end
#Destroy 2 waveguide photon
function waveguide_mul!(result,a::WaveguideDestroy{B,B,2},b,alpha,beta) where {B}
    rmul!(result,beta)
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    add_zerophoton_onephoton!(result,b,alpha*a.factor,timeindex)
    twophotonview = TwophotonTimestepView(b,timeindex,nsteps)
    add_onephoton_twophoton!(result,twophotonview,alpha*a.factor,nsteps)
    return
end


#Create 1 waveguide photon 
function waveguide_mul!(result,a::WaveguideCreate{B,B,1},b,alpha,beta) where {B}
    rmul!(result,beta)
    add_onephoton_zerophoton!(result,b,alpha*a.factor,a.timeindex)
    return
end
#Create 2 waveguide photon
function waveguide_mul!(result,a::WaveguideCreate{B,B,2},b,alpha,beta) where {B}
    rmul!(result,beta)
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    add_onephoton_zerophoton!(result,b,alpha*a.factor,timeindex)
    view = TwophotonTimestepView(result,timeindex,nsteps)
    add_twophoton_onephoton!(view,b,alpha*a.factor)
    return
end

function add_zerophoton_onephoton!(a,b,alpha,timeindex::Int)
    a[1] += alpha*b[timeindex+1]
end

function add_onephoton_zerophoton!(a,b,alpha,timeindex::Int)
    a[1+timeindex] += alpha*b[1]
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
