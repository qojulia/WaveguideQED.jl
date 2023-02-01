
#Abstract type for WaveguideOperators. Used to dispatch special mul! function
abstract type WaveguideOperator <: AbstractOperator{Basis,Basis} end

#Type for annihilation waveguide operator
mutable struct WaveguideDestroy{T} <: WaveguideOperator
    basis_l::Basis
    basis_r::Basis
    factor::ComplexF64
end

#Type for creation annihilation operator
mutable struct WaveguideCreate{T} <:WaveguideOperator
    basis_l::Basis
    basis_r::Basis
    factor::ComplexF64
end

#Methods for copying waveguide operators
function Base.:copy(x::WaveguideDestroy{1})
    WaveguideDestroy{1}(x.basis_l,x.basis_r,x.factor)
end

function Base.:copy(x::WaveguideDestroy{2})
    WaveguideDestroy{2}(x.basis_l,x.basis_r,x.factor)
end

function Base.:copy(x::WaveguideCreate{1})
    WaveguideCreate{1}(x.basis_l,x.basis_r,x.factor)
end

function Base.:copy(x::WaveguideCreate{2})
    WaveguideCreate{2}(x.basis_l,x.basis_r,x.factor)
end


#Method for multiplying, which updates factor in the operator.
function Base.:*(a::Number,b::WaveguideOperator)
    out = copy(b)
    out.factor=out.factor*a
end

Base.:*(b::WaveguideOperator,a::Number)=*(a,b)

function Base.:eltype(x::WaveguideOperator) typeof(x.factor) end


#Methods for tensorproducts between QuantumOptics.jl operators. This is done by forming a LazyTensor.
#TODO: Update method to allow for three or more hilbert spaces.
function QuantumOpticsBase.:tensor(op1::AbstractOperator,op2::WaveguideOperator) 
    btotal = tensor(op1.basis_l,op2.basis_r)
    LazyTensor(btotal,btotal,[1,2],(op1,op2))
end

function QuantumOpticsBase.:tensor(op1::WaveguideOperator,op2::AbstractOperator) 
    btotal = tensor(op1.basis_l,op2.basis_l)
    LazyTensor(btotal,btotal,[1,2],(op1,op2))
end

#Constructors from waveguidebasis.
function QuantumOpticsBase.:destroy(basis::WaveguideBasis{1})
    return WaveguideDestroy{1}(basis,basis,1)
end

function QuantumOpticsBase.:destroy(basis::WaveguideBasis{2})
    return WaveguideDestroy{2}(basis,basis,1)
end

function QuantumOpticsBase.:create(basis::WaveguideBasis{1})
    return WaveguideCreate{1}(basis,basis,1)
end
function QuantumOpticsBase.:create(basis::WaveguideBasis{2})
    return WaveguideCreate{2}(basis,basis,1)
end

function QuantumOpticsBase.:dagger(op::WaveguideCreate)
    @assert op.basis_l == op.basis_r
    destroy(op.basis_l) 
end

function QuantumOpticsBase.:dagger(op::WaveguideDestroy)
    @assert op.basis_l == op.basis_r
    create(op.basis_l) 
end


#Mul! function for CavityWaveguide LazyTensor. 
#TODO: Makes a call to _tp_sum_matmul! which leads to huge allocation. Can be solved by:
#https://github.com/qojulia/QuantumOptics.jl/issues/333
function QuantumOpticsBase.:mul!(result::Ket{B1}, a::LazyTensor{B1,B2,F,I,T}, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis, F,I,T<:Tuple{Vararg{AbstractOperator}}}
    b_data = Base.ReshapedArray(b.data, QuantumOpticsBase._comp_size(basis(b)), ())
    result_data = Base.ReshapedArray(result.data, QuantumOpticsBase._comp_size(basis(result)), ())

    tp_ops = (tuple(( (isa(op,DataOperator) ? op.data : op) for op in a.operators)...), a.indices)
    iso_ops = QuantumOpticsBase._explicit_isometries(a.indices, a.basis_l, a.basis_r)

    QuantumOpticsBase._tp_sum_matmul!(result_data, tp_ops, iso_ops, b_data, alpha * a.factor, beta)
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

function QuantumOpticsBase._tp_sum_get_tmp(op::WaveguideOperator, loc::Integer, arr::AbstractArray{S,N}, sym) where {S,N}
    shp = ntuple(i -> i == loc ? size(op,1) : size(arr,i), N)
    QuantumOpticsBase._tp_matmul_get_tmp(S, shp, sym)
end

function apply_last_op!(result,a,br,α,β)
    for i in axes(br,1)
        waveguide_mul!(view(result,i,:), a, view(br,i,:), α, β)
    end
end

function apply_first_op!(result,a,br,α,β)
    for i in axes(br,2)
        waveguide_mul!(view(result,:,i), a, view(br,:,i), α, β)
    end
end

#Destroy 1 waveguide photon
function waveguide_mul!(result,a::WaveguideDestroy{1},b,alpha,beta)
    rmul!(result,beta)
    add_zerophoton_onephoton!(result,b,alpha,a.basis_l.timeindex)
    return
end

#Destroy 2 waveguide photon
function waveguide_mul!(result,a::WaveguideDestroy{2},b,alpha,beta)
    rmul!(result,beta)
    timeindex = a.basis_l.timeindex
    nsteps = a.basis_l.nsteps
    add_zerophoton_onephoton!(result,b,alpha,timeindex)
    twophotonview = TwophotonView(b,timeindex,nsteps)
    add_onephoton_twophoton!(result,twophotonview,alpha,nsteps)
    return
end


#Create 1 waveguide photon 
function waveguide_mul!(result,a::WaveguideCreate{1},b,alpha,beta)
    rmul!(result,beta)
    add_onephoton_zerophoton!(result,b,alpha,a.basis_l.timeindex)
    return
end

#Create 2 waveguide photon
function waveguide_mul!(result,a::WaveguideCreate{2},b,alpha,beta)
    rmul!(result,beta)
    timeindex = a.basis_l.timeindex
    nsteps = a.basis_l.nsteps
    add_onephoton_zerophoton!(result,b,alpha,timeindex)
    view = TwophotonView(result,timeindex,nsteps)
    add_twophoton_onephoton!(view,b,alpha)
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