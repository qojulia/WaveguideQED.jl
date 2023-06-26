
function Base.:*(a::LazyProduct{B1,B2}, b::AbstractOperator{B1,B2}) where {B1,B2}
    LazyProduct((a.operators..., b), a.factor)
end

function Base.:*(a::AbstractOperator{B1,B2}, b::LazyProduct{B1,B2}) where {B1,B2}
    LazyProduct((a, b.operators...), b.factor)
end

function Base.:*(a::LazyProduct{B1,B2}, b::Operator{B1,B2}) where {B1,B2}
    LazyProduct((a.operators..., b), a.factor)
end

function Base.:*(a::Operator{B1,B2}, b::LazyProduct{B1,B2}) where {B1,B2}
    LazyProduct((a, b.operators...), b.factor)
end

function Base.:*(a::LazySum{B1,B2}, b::AbstractOperator{B1,B2}) where {B1,B2}
    QuantumOpticsBase.check_samebases(a,b)
    LazySum(a.basis_l, a.basis_r, factors, a.operators .* b)
end

function Base.:*(a::LazySum{B1,B2}, b::LazySum{B1,B2}) where {B1,B2}
    QuantumOpticsBase.check_samebases(a,b)
    c = Array{AbstractOperator}(undef,length(a.operators)*length(b.operators))
    factors = similar(a.factors,length(a.operators)*length(b.operators))
    k = 1
    for i in eachindex(a.operators)
        for j in eachindex(b.operators)
            factors[k] = a.factors[i] * b.factors[j]
            c[k] = a.operators[i] * b.operators[j]
            k += 1
        end
    end
    @samebases LazySum(a.basis_l, a.basis_r, factors, (c...,))
end