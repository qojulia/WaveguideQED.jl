abstract type  CavityWaveguideOperator{BL,BR} <: AbstractOperator{BL,BR} end


#CavityWaveguide operator for performance enhancement
mutable struct CavityWaveguideAbsorption{B1,B2} <: CavityWaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    op::AbstractOperator
    loc
end

#Type for creation annihilation operator
mutable struct CavityWaveguideEmission{B1,B2} <: CavityWaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    op::AbstractOperator
    loc
end

function Base.:copy(x::CavityWaveguideEmission)
    CavityWaveguideEmission(x.basis_l,x.basis_r,x.factor,x.op,x.loc)
end

function Base.:copy(x::CavityWaveguideAbsorption)
    CavityWaveguideAbsorption(x.basis_l,x.basis_r,x.factor,x.op,x.loc)
end

#Method for multiplying, which updates factor in the operator.
function Base.:*(a::Number,b::CavityWaveguideOperator)
    out = copy(b)
    out.factor=out.factor*a
end

Base.:*(b::CavityWaveguideOperator,a::Number)=*(a,b)

function Base.:eltype(x::CavityWaveguideOperator) typeof(x.factor) end

function get_cavity_operator(a::CavityWaveguideEmission)
    destroy(a.basis_l.bases[a.loc[2]])
end

function get_cavity_operator(a::CavityWaveguideAbsorption)
    create(a.basis_l.bases[a.loc[2]])
end


#Constructors from a Cavity and Waveguide basis
function absorption(b1::WaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),destroy(b1),[1,2])
end

function emission(b1::WaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),create(b1),[1,2])
end

function absorption(b1::FockBasis,b2::WaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),destroy(b2),[2,1])
end

function emission(b1::FockBasis,b2::WaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),create(b2),[2,1])
end


function QuantumOpticsBase.:tensor(op1::AbstractOperator,op2::T) where {T<:CavityWaveguideAbsorption}
    btotal = tensor(op1.basis_l,op2.basis_r)
    if isequal(op1,identityoperator(basis(op1)))
        CavityWaveguideAbsorption(btotal,btotal,op2.factor,op2.op,op2.loc .+length(basis(op1).shape))
    else
        sorted_idx = sortperm([1,op2.loc[1]+1,op2.loc[2]+1])
        LazyTensor(btotal,btotal,[1,op2.loc[1]+1,op2.loc[2]+1][sorted_idx],(op1,op2.op,get_cavity_operator(op2))[sorted_idx])
    end
end

function QuantumOpticsBase.:tensor(op1::T,op2::AbstractOperator) where {T<:CavityWaveguideAbsorption}
    btotal = tensor(op1.basis_l,op2.basis_r)
    if isequal(op2,identityoperator(basis(op2)))
        CavityWaveguideAbsorption(btotal,btotal,op1.factor,op1.op,op1.loc)
    else
        sorted_idx = sortperm([op1.loc[1]+1,op1.loc[2]+1,length(btotal.shape)])
        LazyTensor(btotal,btotal,[op1.loc[1]+1,op1.loc[2]+1,length(btotal.shape)][sorted_idx],(op1.op,get_cavity_operator(op1),op2)[sorted_idx])
    end
end

function QuantumOpticsBase.:tensor(op1::AbstractOperator,op2::T) where {T<:CavityWaveguideEmission}
    btotal = tensor(op1.basis_l,op2.basis_r)
    if isequal(op1,identityoperator(basis(op1)))
        CavityWaveguideEmission(btotal,btotal,op2.factor,op2.op,op2.loc .+length(basis(op1).shape))
    else
        sorted_idx = sortperm([1,op2.loc[1]+1,op2.loc[2]+1])
        LazyTensor(btotal,btotal,[1,op2.loc[1]+1,op2.loc[2]+1][sorted_idx],(op1,op2.op,get_cavity_operator(op2))[sorted_idx])
    end
end

function QuantumOpticsBase.:tensor(op1::T,op2::AbstractOperator) where {T<:CavityWaveguideEmission}
    btotal = tensor(op1.basis_l,op2.basis_r)
    if isequal(op2,identityoperator(basis(op2)))
        CavityWaveguideEmission(btotal,btotal,op1.factor,op1.op,op1.loc)
    else
        sorted_idx = sortperm([op1.loc[1]+1,op1.loc[2]+1,length(btotal.shape)])
        LazyTensor(btotal,btotal,[op1.loc[1]+1,op1.loc[2]+1,length(btotal.shape)][sorted_idx],(op1.op,get_cavity_operator(op1),op2)[sorted_idx])
    end
end


function QuantumOpticsBase.:mul!(result::Ket{B1}, a::CavityWaveguideOperator, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis}
    if a.loc[1] == 1
        return matmul_first!(result.data, a, b.data, alpha, beta)
    elseif a.loc[1] == length(b.basis.shape)
        return matmul_last!(result.data, a, b.data, alpha, beta)
    end
    error("Waveguide operators in CavityWaveguide operators has to be first or last")
end


function QuantumOpticsBase.:mul!(result::Bra{B1}, a::Bra{B2}, b::CavityWaveguideOperator, alpha, beta) where {B1<:Basis,B2<:Basis}
    if b.loc[1] == 1
        return matmul_first!(result.data, b, a.data, alpha, beta)
    elseif b.loc[1] == length(a.basis.shape)
        return matmul_last!(result.data, b, a.data, alpha, beta)
    end
    error("Waveguide operators in CavityWaveguide operators has to be first or last")
end

   
#Called from _tp_matmul!
#Calls CavityWaveguide operator  on correct view of subsets of the state.
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

#Same as _tp_matmul_first! But indexed in another way.
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

function waveguide_mul_last!(result, a::CavityWaveguideAbsorption, b, alpha::Number, beta::Number)
    b_data = reshape(b, :, length(basis(a).bases[end]))
    result_data = reshape(result,  :, length(basis(a).bases[end]))
    for i in 2:size(b_data,1)
        waveguide_mul!(view(result_data,i,:),a.op,view(b_data,i-1,:),sqrt(i-1)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,1,:),beta)
end


function waveguide_mul_first!(result, a::CavityWaveguideAbsorption, b, alpha::Number, beta::Number)
    b_data = reshape(b, length(basis(a).bases[1]),:)
    result_data = reshape(result,length(basis(a).bases[1]),:)
    for i in 2:size(b_data,2)
        waveguide_mul!(view(result_data,:,i),a.op,view(b_data,:,i-1),sqrt(i-1)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,:,1),beta)
end


function waveguide_mul_last!(result, a::CavityWaveguideEmission, b, alpha::Number, beta::Number)
    b_data = reshape(b, :, length(basis(a).bases[end]))
    result_data = reshape(result,  :, length(basis(a).bases[end]))
    for i in 1:size(b_data,1)-1
        waveguide_mul!(view(result_data,i,:),a.op,view(b_data,i+1,:),sqrt(i)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,size(b_data,1),:),beta)
end


function waveguide_mul_first!(result, a::CavityWaveguideEmission, b, alpha::Number, beta::Number)
    b_data = reshape(b,length(basis(a).bases[1]),:)
    result_data = reshape(result,length(basis(a).bases[1]),:)
    for i in 1:size(b_data,2)-1
        waveguide_mul!(view(result_data,:,i),a.op,view(b_data,:,i+1),sqrt(i)*alpha*a.factor,beta)
    end
    rmul!(view(result_data,:,size(b_data,2)),beta)
end


function apply_last_op!(result,a::CavityWaveguideOperator,br,α,β)
    for i in axes(br,1)
        waveguide_mul_last!(view(result,i,:), a, view(br,i,:), α, β)
    end
end

function apply_first_op!(result,a::CavityWaveguideOperator,br,α,β)
    for i in axes(br,2)
        waveguide_mul_first!(view(result,:,i), a, view(br,:,i), α, β)
    end
end