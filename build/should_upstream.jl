
function QuantumOpticsBase.:+(a::LazyTensor{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    LazySum(a,b)
end

function QuantumOpticsBase.:-(a::LazyTensor{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    LazySum([1,-1],[a,b])
end

function QuantumOpticsBase.:+(a::CavityWaveguideOperator,b::CavityWaveguideOperator)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum(a,b)
end

function QuantumOpticsBase.:-(a::CavityWaveguideOperator,b::CavityWaveguideOperator)
    @assert a.basis_l == b.basis_l
    @assert a.basis_r == b.basis_r
    LazySum([1,-1],[a,b])
end

function QuantumOpticsBase.:+(a::LazyTensor{B1,B2},b::Operator{B1,B2}) where {B1,B2}
    LazySum(a) + b
end
function QuantumOpticsBase.:+(a::Operator{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    +(b,a)
end

function QuantumOpticsBase.:⊗(a::LazyTensor,b::Operator)
    a ⊗ LazyTensor(b.basis_l,b.basis_r,[1],(b,),1)
end

function QuantumOpticsBase.:⊗(a::Operator,b::LazyTensor)
    LazyTensor(a.basis_l,a.basis_r,[1],(a,),1) ⊗ b
end

function QuantumOpticsBase.:⊗(a::Operator,b::LazySum)
    btotal_l = a.basis_l ⊗ b.basis_l
    btotal_r = a.basis_r ⊗ b.basis_r
    ops = Vector{AbstractOperator}(undef,length(b.operators))
    for i in eachindex(ops)
        ops[i] = a ⊗ b.operators[i]
    end
    LazySum(btotal_l,btotal_r,b.factors,ops)
end
function QuantumOpticsBase.:⊗(a::LazySum,b::Operator)
    btotal_l = a.basis_l ⊗ b.basis_l
    btotal_r = a.basis_r ⊗ b.basis_r
    ops = Vector{AbstractOperator}(undef,length(a.operators))
    for i in eachindex(ops)
        ops[i] = a.operators[i] ⊗ b
    end
    LazySum(btotal_l,btotal_r,a.factors,ops)
end


function QuantumOpticsBase.:⊗(a::Operator,b::LazyProduct)
    btotal_l = a.basis_l ⊗ b.basis_l
    btotal_r = a.basis_r ⊗ b.basis_r
    ops = Vector{AbstractOperator}(undef,length(b.operators))
    for i in eachindex(ops)
        ops[i] = a ⊗ b.operators[i]
    end
    LazyProduct(btotal_l,btotal_r,b.factors,ops)
end
function QuantumOpticsBase.:⊗(a::LazyProduct,b::Operator)
    btotal_l = a.basis_l ⊗ b.basis_l
    btotal_r = a.basis_r ⊗ b.basis_r
    ops = Vector{AbstractOperator}(undef,length(a.operators))
    for i in eachindex(ops)
        ops[i] = a.operators[i] ⊗ b
    end
    LazyProduct(btotal_l,btotal_r,a.factors,ops)
end

function QuantumOpticsBase.:⊗(a::AbstractOperator,b::LazyTensor)
    btotal_l = a.basis_l ⊗ b.basis_l
    btotal_r = a.basis_r ⊗ b.basis_r
    indices = (1,(b.indices .+ 1)...)
    ops = (a,b.operators...)
    LazyTensor(btotal_l,btotal_r,indices,ops,b.factor)
end
function QuantumOpticsBase.:⊗(a::LazyTensor,b::AbstractOperator)
    btotal_l = a.basis_l ⊗ b.basis_l
    btotal_r = a.basis_r ⊗ b.basis_r
    indices = (a.indices...,length(a.indices)+1)
    ops = (a.operators...,b)
    LazyTensor(btotal_l,btotal_r,indices,ops,a.factor)
end
#Commented out as this needs furtherwork
"""

function Base.:+(a::LazySum{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    LazySum([a.factors; 1], (a.operators..., b))    
end

function Base.:+(a::LazyTensor{B1,B2},b::LazySum{B1,B2}) where {B1,B2}
    +(b,a)   
end

function Base.:-(a::LazySum{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    LazySum([a.factors; -1], (a.operators..., b))   
end

function Base.:-(a::LazyTensor{B1,B2},b::LazySum{B1,B2}) where {B1,B2}
    LazySum([1; b.factors...], (a, b.operators...))
end

function QuantumOpticsBase.:+(a::LazySum{B1,B2},b::Operator{B1,B2}) where {B1,B2}
    LazySum([a.factors; 1], (a.operators..., b))      
end

QuantumOpticsBase.:+(a::Operator,b::LazySum) = +(b,a)

function QuantumOpticsBase.:-(a::LazySum{B1,B2},b::Operator{B1,B2}) where {B1,B2}
    LazySum([a.factors; -1], (a.operators..., b))       
end

function QuantumOpticsBase.:-(a::Operator,b::LazySum{B1,B2}) where {B1,B2}
    LazySum([1; -b.factors...], (a, b.operators...))       
end
"""