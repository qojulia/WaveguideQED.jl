
function QuantumOpticsBase.:+(a::LazyTensor{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    LazySum(a,b)
end

function QuantumOpticsBase.:-(a::LazyTensor{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    LazySum([1,-1],[a,b])
end

function QuantumOpticsBase.:+(a::LazyTensor{B1,B2},b::Operator{B1,B2}) where {B1,B2}
    LazySum(a) + b
end
function QuantumOpticsBase.:+(a::Operator{B1,B2},b::LazyTensor{B1,B2}) where {B1,B2}
    +(b,a)
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
    LazySum([1; b.factors...], (a, b.operaterors...))
end

function QuantumOpticsBase.:+(a::LazySum{B1,B2},b::Operator{B1,B2}) where {B1,B2}
    LazySum([a.factors; 1], (a.operators..., b))      
end

QuantumOpticsBase.:+(a::Operator,b::LazySum) = +(b,a)

function QuantumOpticsBase.:-(a::LazySum{B1,B2},b::Operator{B1,B2}) where {B1,B2}
    LazySum([a.factors; -1], (a.operators..., b))       
end

function QuantumOpticsBase.:-(a::Operator,b::LazySum{B1,B2}) where {B1,B2}
    LazySum([1; -b.factors...], (a, b.operaterors...))       
end
"""