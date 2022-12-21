#Commented out as this needs furtherwork
"""
function Base.:+(op1::LazyTensor{B1,B2},op2::LazyTensor{B1,B2}) where {B1,B2}
    LazySum(op1,op2)
end


function Base.:-(op1::LazyTensor{B1,B2},op2::LazyTensor{B1,B2}) where {B1,B2}
    LazySum([1,-1],(op1,op2))
end

function Base.:+(op1::LazySum{B1,B2},op2::LazyTensor{B1,B2}) where {B1,B2}
    op1+LazySum(op2)    
end

function Base.:+(op1::LazyTensor{B1,B2},op2::LazySum{B1,B2}) where {B1,B2}
    +(op2,op1)   
end

function Base.:-(op1::LazySum{B1,B2},op2::LazyTensor{B1,B2}) where {B1,B2}
    op1-LazySum(op2)   
end

function Base.:-(op1::LazyTensor{B1,B2},op2::LazySum{B1,B2}) where {B1,B2}
    LazySum(op1)-op2
end

function Base.:+(op1::LazySum{B1,B2},op2::Operator) where {B1,B2}
    op1+LazySum(op2)      
end

Base.:+(op1::Operator,op2::LazySum) = +(op2,op1)

function Base.:-(op1::LazySum{B1,B2},op2::Operator) where {B1,B2}
    op1-LazySum(op2)       
end

function Base.:-(op1::Operator,op2::LazySum{B1,B2}) where {B1,B2}
    LazySum(op1)-op2       
end

"""