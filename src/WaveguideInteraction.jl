mutable struct WaveguideInteraction{B1,B2} <: AbstractOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    op1::AbstractOperator
    op2::AbstractOperator
    loc
end

function Base.:eltype(x::WaveguideInteraction) typeof(x.factor) end
function Base.:copy(x::WaveguideInteraction)
    WaveguideInteraction(x.basis_l,x.basis_r,x.factor,x.op1,x.op2,x.loc)
end

Base.:*(x::WaveguideInteraction{BL,BR},y::WaveguideInteraction{BL,BR}) where {BL,BR} = LazyProduct((x,y),x.factor*y.factor)

identityoperator(x::WaveguideInteraction) = identityoperator(x.basis_l)
function identityoperator(::Type{T},::Type{ComplexF64}, b1::Basis, b2::Basis) where {T<:WaveguideInteraction}
    @assert b1==b2
    identityoperator(b1)
end
#Method for multiplying, which updates factor in the operator.
function Base.:*(a::Number,b::WaveguideInteraction)
    out = copy(b)
    out.factor=out.factor*a
    out
end
Base.:*(b::WaveguideInteraction,a::Number)=*(a,b)

function QuantumOpticsBase.:+(a::WaveguideInteraction,b::WaveguideInteraction)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum(a,b)
end
function QuantumOpticsBase.:-(a::WaveguideInteraction,b::WaveguideInteraction)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum([1,-1],[a,b])
end

function set_waveguidetimeindex!(op::WaveguideInteraction,timeindex)
    op.op1.timeindex = timeindex
    op.op2.timeindex = timeindex
end
function get_waveguide_operators(op::WaveguideInteraction)
    [op.op1,op.op2]
end

function tensor(a::AbstractOperator,b::WaveguideInteraction)
    btotal = tensor(a.basis_l,b.basis_r)
    if isequal(a,identityoperator(basis(a)))
        WaveguideInteraction(btotal,btotal,b.factor,b.op1,b.op2,b.loc .+length(basis(a).shape))
    else
        sorted_idx = sortperm([1,b.loc[1]+1,b.loc[2]+1])
        LazyTensor(btotal,btotal,[1,b.loc[1]+1,b.loc[2]+1][sorted_idx],(a,b)[sorted_idx])
    end
end
function tensor(a::WaveguideInteraction,b::AbstractOperator)
    btotal = tensor(a.basis_l,b.basis_r)
    if isequal(b,identityoperator(basis(b)))
        WaveguideInteraction(btotal,btotal,a.factor,a.op1,a.op2,a.loc)
    else
        sorted_idx = sortperm([a.loc[1]+1,a.loc[2]+1,length(btotal.shape)])
        LazyTensor(btotal,btotal,[a.loc[1]+1,a.loc[2]+1,length(btotal.shape)][sorted_idx],(a,b)[sorted_idx])
    end
end

function mul!(result::Ket{B1}, a::WaveguideInteraction{B1,B2}, b::Ket{B2}, alpha, beta) where {B1<:Basis,B2<:Basis}
    dims = basis(result).shape
    
    alldims = [1:length(dims)...]
    i = a.loc[1]
    exclude_dims = [i]
    otherdims = setdiff(alldims, exclude_dims)
    iter = (dims[otherdims]...,)

    idx_vec1=ones(Int64,length(dims))
    idx_vec2=ones(Int64,length(dims))
    
    range_idx = Vector{Int}(undef, length(iter))
    range_idx .= otherdims

    idx_vec1[i] = 1
    idx_vec2[i] = 2
    stride = linear_index(dims,idx_vec2)- linear_index(dims,idx_vec1)
    idx_vec2[i] = 1
    if length(dims) == 1
        li1 = linear_index(dims,idx_vec1)
        end_idx = dims[i]
        waveguide_interaction_mul!(view(result.data,li1:stride:li1+(end_idx-1)*stride),a.op1,a.op2,view(b.data,li1:stride:li1+(end_idx-1)*stride),alpha,beta)            
    else
        iterate_over_interaction!(result.data,b.data,a,a.factor*alpha,beta,iter,idx_vec1,range_idx,1,stride,dims[i],dims)
    end

    return result
end

function iterate_over_interaction!(result,b,a,alpha,beta,iter::Tuple,idx_vec1,range_idx,idx,stride,end_idx,dims)
    if length(iter) == idx
        for i in Base.OneTo(iter[idx])
            @inbounds idx_vec1[range_idx[idx]] = i
            li1 = linear_index(dims,idx_vec1)
            @uviews result b begin
            waveguide_interaction_mul!(view(result,li1:stride:li1+(end_idx-1)*stride),a.op1,a.op2,view(b,li1:stride:li1+(end_idx-1)*stride),alpha,beta)            
            end
        end
    else
        for j in Base.OneTo(iter[idx])
            @inbounds idx_vec1[range_idx[idx]] = j
            iterate_over_interaction!(result,b,a,alpha,beta,iter,idx_vec1,range_idx,idx+1,stride,end_idx,dims)
        end
    end
end

function waveguide_interaction_mul!(result,a::WaveguideCreate{B1,B2,1,idx1},bop::WaveguideDestroy{B1,B2,1,idx2},b,alpha,beta) where {B1,B2,idx1,idx2}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    
    @inbounds result[1+(idx1-1)*nsteps+timeindex] += alpha*a.factor*bop.factor*b[1+(idx2-1)*nsteps+timeindex] 
    result
end

function waveguide_interaction_mul!(result,a::WaveguideCreate{B1,B2,2,idx1},bop::WaveguideDestroy{B1,B2,2,idx2},b,alpha,beta) where {B1,B2,idx1,idx2}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw = get_number_of_waveguides(a.basis_l)    

    @inbounds result[1+(idx1-1)*nsteps+timeindex] += alpha*a.factor*bop.factor*b[1+(idx2-1)*nsteps+timeindex] 
    
    twophotonview_b = TwoPhotonTimestepView(b,timeindex,nsteps,Nw*nsteps+1+(idx2-1)*(nsteps*(nsteps+1))÷2)
    i,j = min(idx1,idx2),max(idx1,idx2)
    index = (i-1)*Nw + j - (i*(i+1))÷2
    twophotonview_r = TwoWaveguideTimestepView(result,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx1)
    axpy!(alpha*a.factor*bop.factor,twophotonview_b,twophotonview_r)
    
    i,j = min(idx1,idx2),max(idx1,idx2)
    index = (i-1)*Nw + j - (i*(i+1))÷2
    twophotonview_b = TwoWaveguideTimestepView(b,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx2)
    twophotonview_r = TwoPhotonTimestepView(result,timeindex,nsteps,1+Nw*nsteps+(idx1-1)*(nsteps*(nsteps+1))÷2)
    axpy!(alpha*a.factor*bop.factor,twophotonview_b,twophotonview_r)
    
    @simd for k in filter(x -> x != idx1 &&  x != idx2, 1:Nw)
        m,l = min(k,idx1),max(k,idx1)
        index_r = (m-1)*Nw + l - (m*(m+1))÷2
        twophotonview_r = TwoWaveguideTimestepView(result,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index_r-1)*nsteps^2,m==idx1)
        
        i,j = min(k,idx2),max(k,idx2)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        twophotonview_b = TwoWaveguideTimestepView(b,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx2)
        axpy!(alpha*a.factor*bop.factor,twophotonview_b,twophotonview_r)
    end
    result
end

function waveguide_interaction_mul!(result,a::WaveguideCreate{B1,B2,2,idx},bop::WaveguideDestroy{B1,B2,2,idx},b,alpha,beta) where {B1,B2,idx}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw = get_number_of_waveguides(a.basis_l)

    @inbounds result[1+(idx-1)*nsteps+timeindex] += alpha*a.factor*bop.factor*b[1+(idx-1)*nsteps+timeindex] 
    
    twophotonview_b = TwoPhotonTimestepView(b,timeindex,nsteps,Nw*nsteps+1+(idx-1)*(nsteps*(nsteps+1))÷2)
    twophotonview_r = TwoPhotonTimestepView(result,timeindex,nsteps,Nw*nsteps+1+(idx-1)*(nsteps*(nsteps+1))÷2)
    axpy!(alpha*a.factor*bop.factor,twophotonview_b,twophotonview_r)
    
    @simd for k in filter(x -> x != idx, 1:Nw)
        m,l = min(k,idx),max(k,idx)
        index_r = (m-1)*Nw + l - (m*(m+1))÷2
        twophotonview_r = TwoWaveguideTimestepView(result,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index_r-1)*nsteps^2,m==idx)
        twophotonview_b = TwoWaveguideTimestepView(b,timeindex,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index_r-1)*nsteps^2,m==idx)
        axpy!(alpha*a.factor*bop.factor,twophotonview_b,twophotonview_r)
    end
    result
end

function waveguide_interaction_mul!(result,a::WaveguideDestroy{B1,B2,2,idx},bop::WaveguideDestroy{B1,B2,2,idx},b,alpha,beta) where {B1,B2,idx}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw = get_number_of_waveguides(a.basis_l)

    result[1] += sqrt(2)*alpha*a.factor*bop.factor*b[1+Nw*nsteps+(idx-1)*(nsteps*(nsteps+1))÷2+twophoton_index(timeindex,nsteps,timeindex)] 
    result
end
function waveguide_interaction_mul!(result,a::WaveguideCreate{B1,B2,2,idx},bop::WaveguideCreate{B1,B2,2,idx},b,alpha,beta) where {B1,B2,idx}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw = get_number_of_waveguides(a.basis_l)
    result[1+Nw*nsteps+(idx-1)*(nsteps*(nsteps+1))÷2+twophoton_index(timeindex,nsteps,timeindex)] += sqrt(2)*alpha*a.factor*bop.factor*b[1] 
    result
end

function waveguide_interaction_mul!(result,a::WaveguideDestroy{B1,B2,2,idx1},bop::WaveguideDestroy{B1,B2,2,idx2},b,alpha,beta) where {B1,B2,idx1,idx2}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw = get_number_of_waveguides(a.basis_l)
    m,l = min(idx1,idx2),max(idx1,idx2)
    index_r = (m-1)*Nw + l - (m*(m+1))÷2
    result[1] += alpha*a.factor*bop.factor*b[1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index_r-1)*nsteps^2+(timeindex-1)*nsteps + timeindex] 
    result
end
function waveguide_interaction_mul!(result,a::WaveguideCreate{B1,B2,2,idx1},bop::WaveguideCreate{B1,B2,2,idx2},b,alpha,beta) where {B1,B2,idx1,idx2}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    Nw = get_number_of_waveguides(a.basis_l)
    m,l = min(idx1,idx2),max(idx1,idx2)
    index_r = (m-1)*Nw + l - (m*(m+1))÷2
    result[1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index_r-1)*nsteps^2+(timeindex-1)*nsteps + timeindex]  += alpha*a.factor*bop.factor*b[1] 
end




function waveguide_interaction_mul!(result,a::WaveguideDestroy{B1,B2,Np,idx1},bop::WaveguideCreate{B1,B2,1,idx2},b,alpha,beta) where {B1,B2,Np,idx1,idx2}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    
    @inbounds result[1] = alpha*a.factor*bop.factor*b[1]
    @inbounds result[1+(idx1-1)*nsteps+timeindex] += alpha*a.factor*bop.factor*b[1+(idx2-1)*nsteps+timeindex] 
    result
end


function Base.:*(a::WaveguideCreate{B1,B2,Np,idx1},b::WaveguideDestroy{B1,B2,Np,idx2}) where {B1,B2,Np,idx1,idx2}
    WaveguideInteraction(basis(a),basis(a),complex(1.0),a,b,1)
end
function Base.:*(a::T1,b::T2) where {T1<: WaveguideOperatorT,T2<: WaveguideOperatorT}
    @assert a.indices[1] == b.indices[1]
    WaveguideInteraction(basis(a),basis(a),a.factor*b.factor,a.operators[1],b.operators[1],a.indices[1])
end

const TensorWaveguideInteraction = LazyTensor{B,B,F,V,T} where {B,F,V,T<:Tuple{WaveguideInteraction}}


function expect(a::T,psi::Ket) where T<:Union{WaveguideInteraction,TensorWaveguideInteraction}
    out = 0
    for i in 1:get_nsteps(basis(a))
        #println(i)
        out += expect(a,psi,i)
    end
    out
end

function expect(a::T,psi::Ket,i) where T<:Union{WaveguideInteraction,TensorWaveguideInteraction}
    set_waveguidetimeindex!(a,i)
    dot(psi.data,(a*psi).data)
end