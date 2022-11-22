

import Base:+,-,/,*,zero,copy

export annihilation,creation,adw,wda,number,CWOperator,CWLazySum,CWLazyProduct,CWSingleOperator,AdW,WdA,Creation,Annihilation,ParticleNumber

abstract type CWOperator end
mutable struct AdW <:CWOperator
    factor::ComplexF64
end
mutable struct WdA <:CWOperator 
    factor::ComplexF64
end
mutable struct Creation <:CWOperator 
    factor::ComplexF64
end
mutable struct Annihilation <: CWOperator
    factor::ComplexF64
end
mutable struct ParticleNumber <: CWOperator
    factor::ComplexF64
end
mutable struct Identity <:CWOperator
    factor::ComplexF64
end

"""
mutable struct CWSingleOperator <: CWOperator
    factor::ComplexF64
    operation::Function
end
"""

copy(operator::T) where T<:CWOperator = T(operator.factor)

function *(operator::CWOperator,a::ComplexF64)
    cpy = copy(operator)
    cpy.factor = cpy.factor*a
    return cpy
end
*(operator::CWOperator,a::Number) = *(operator,complex(float(a)))
*(a::Number,operator::CWOperator) = *(operator,a)
/(operator::CWOperator,a::Number) = *(operator,complex(float(1/a)))

mul!(state::CWState,operator::AdW) = adw!(state;g=operator.factor) 
mul!(state::CWState,operator::WdA) = wda!(state;g=operator.factor) 
mul!(state::CWState,operator::Creation) = creation!(state;g=operator.factor) 
mul!(state::CWState,operator::Annihilation) = annihilation!(state;g=operator.factor) 
mul!(state::CWState,operator::ParticleNumber) = number!(state;g=operator.factor) 
#mul!(state::CWState,operator::Identiy) = identity!(state;g=operator.factor)
mul!(state::CWState,operator::Identity) = "Do nothing"

"""
function *(operator::CWSingleOperator,state::CWState) 
    operator.operation(state,g=operator.factor)
end
"""

mutable struct CWLazySum <:CWOperator
    factor::ComplexF64
    operators::Array{CWOperator}
end

copy(operator::CWLazySum) = CWLazySum(operator.factor,copy(operator.operators))

function copy(vec::Vector{CWOperator})
    out = Array{CWOperator}(undef,length(vec))
    for (i,O) in enumerate(vec)
        out[i] = O
    end
    return out
end

+(operator1::CWOperator,operator2::CWOperator) = CWLazySum(complex(1.0),[operator1,operator2])

function +(operator1::CWLazySum,operator2::CWOperator)
    cpy = copy(operator1)
    for i in 1:length(cpy.operators)
        cpy.operators[i] = cpy.factor*cpy.operators[i] 
    end
    cpy.factor = 1
    push!(cpy.operators,operator2)
    return cpy
end
+(operator1::CWOperator,operator2::CWLazySum) = +(operator2,operator1)

function +(operator1::CWLazySum,operator2::CWLazySum)
    cpy1 = copy(operator1)
    for i in 1:length(operator2.operators) 
        cpy1.operators[i] = cpy1.factor*cpy1.operators[i] 
        push!(cpy1.operators,operator2.factor*operator2.operators[i])
    end
    cpy1.factor = 1
    return cpy1
end

function -(operator1::CWOperator,operator2::CWOperator) 
    cpy = copy(operator2)
    cpy.factor = -cpy.factor
    CWLazySum(complex(1.0),[operator1,cpy])
end

function -(operator1::CWLazySum,operator2::CWOperator) 
    cpy = copy(operator2)
    cpy.factor = -cpy.factor
    +(operator1,operator2)
end

function mul!(state::CWState,operator::CWLazySum)
    for O in operator.operators
        mul!(state,operator.factor*O)
    end
end


mutable struct CWLazyProduct <: CWOperator
    factor::ComplexF64
    operators::Array{CWOperator}
end

copy(operator::CWLazyProduct) = CWLazyProduct(operator.factor,copy(operator.operators))


function *(operator1::CWOperator,operator2::CWOperator)
    cpy1 = copy(operator1)
    cpy2 = copy(operator2)
    cpy1.factor = 1
    cpy2.factor = 1
    CWLazyProduct(operator1.factor*operator1.factor,[cpy1,cpy2])
end

function *(operator1::CWLazyProduct,operator2::CWOperator)
    cpy1 = copy(operator1)
    cpy1 = operator2.factor*cpy1
    cpy2 = copy(operator2)
    cpy2.factor = 1
    push!(cpy1.operators,cpy2)
    return cpy1
end

function *(operator1::CWOperator,operator2::CWLazyProduct)
    cpy2 = copy(operator2)
    cpy2 = operator2.factor*cpy2
    cpy1 = copy(operator1)
    cpy1.factor = 1
    pushfirst!(cpy2.operators,cpy1)
    return cpy2
end

function *(operator1::CWLazyProduct,operator2::CWLazySum)
    cpy = copy(operator2)
    cpy = operator1.factor*cpy
    for (i,O_sum) in enumerate(cpy.operators)
        for O_prod in reverse(operator1.operators)
            O_sum = O_prod*O_sum
        end
        cpy.operators[i] = O_sum
    end
    return cpy
end

function *(operator1::CWLazySum,operator2::CWLazyProduct)
    cpy = copy(operator1)
    cpy = operator2.factor*cpy
    for (i,O_sum) in enumerate(cpy.operators)
        for O_prod in operator2.operators
            O_sum = O_sum*O_prod
        end
        cpy.operators[i] = O_sum
    end
    return cpy
end


function *(operator1::CWLazyProduct,operator2::CWLazyProduct)
    cpy = copy(operator1)
    cpy = operator2.factor*cpy
    for O in operator2.operators
        push!(cpy.operators,O)
    end
    return cpy
end

function mul!(state::CWState,operator::CWLazyProduct)
    save_state!(state)
    save_diff!(state)
    zero_diff!(state)
    O = operator.factor*operator.operators[end]
    mul!(state,O)
    set_equal_diff!(state)
    for O in operator.operators[end-1:-1:2]
        mul!(state,O)
        set_equal_diff!(state)
    end
    mul!(state,operator.operators[1])
    load_diff!(state)
    load_state!(state)
end



mutable struct CWLazyTensor <: CWOperator
    cw_operator::CWOperator
    qo_opeartor::Operator
end

function lazytensor(op1::T,op2::Operator) where T<:CWOperator
    return CWLazyTensor(op1,op2)
end

function lazytensor(op1::CWLazyProduct,op2::Operator)
    cpy = copy(op1)
    for (i,O) in enumerate(cpy.operators)
        cpy.operators[i] = lazytensor(O,op2)
    end
    return cpy
end

function lazytensor(op1::CWLazySum,op2::Operator)
    cpy = copy(op1)
    for (i,O) in enumerate(cpy.operators)
        cpy.operators[i] = lazytensor(O,op2)
    end
    return cpy
end

⊗(a::CWOperator,b::Operator) = lazytensor(a,b)
⊗(a::Operator,b::CWOperator) = lazytensor(b,a)


function mul!(state::CWTensorState,operator::CWLazyTensor)
    state.ketbuffer = operator.qo_opeartor*state.ket
    for i in 1:length(state.cwstates)
        mul!(state.cwstates[i],operator.cw_operator)
    end
end

function annihilation(b::cwbasis)
    return Annihilation(1)
end
function creation(b::cwbasis)
    return Creation(1)
end
function adw(b::cwbasis)
    return AdW(1)
end
function wda(b::cwbasis)
    return WdA(1)
end
function number(b::cwbasis)
    return ParticleNumber(1)
end

annihilation(state::CWState) = annihilation(state.basis)
creation(state::CWState) = creation(state.basis)
adw(state::CWState) = adw(state.basis)
wda(state::CWState) = wda(state.basis)
number(state::CWState) = number(state.basis)

