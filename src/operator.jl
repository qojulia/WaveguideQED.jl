
import Base:+,-,/,*,zero

export annihilation,creation,adw,wda,number,CWOperator,CWLazySum,CWLazyProduct,CWSingleOperator

abstract type CWOperator end

mutable struct CWSingleOperator <: CWOperator
    factor::ComplexF64
    operation::Function
end

copy(operator::CWSingleOperator) = CWSingleOperator(operator.factor,operator.operation)

function *(operator::CWOperator,a::ComplexF64)
    cpy = copy(operator)
    cpy.factor = cpy.factor*a
    return cpy
end
*(operator::CWOperator,a::Number) = *(operator,complex(float(a)))
*(a::Number,operator::CWOperator) = *(operator,a)
/(operator::CWOperator,a::Number) = *(operator,complex(float(1/a)))


function *(operator::CWSingleOperator,state::cwstate) 
    operator.operation(state,g=operator.factor)
end


mutable struct CWLazySum <:CWOperator
    factor::ComplexF64
    operators::Array{CWOperator}
end

copy(operator::CWLazySum) = CWLazySum(operator.factor,operator.operators)


+(operator1::CWOperator,operator2::CWOperator) = CWLazySum(complex(1.0),[operator1 operator2])
function -(operator1::CWOperator,operator2::CWOperator) 
    cpy = copy(operator2)
    cpy.factor = -cpy.factor
    CWLazySum(complex(1.0),[operator1 cpy])
end

function *(operator::CWLazySum,state::cwstate)
    for O in operator.operators
        O.factor = O.factor*operator.factor 
        O*state
        O.factor = O.factor/operator.factor
    end
end

mutable struct CWLazyProduct <: CWOperator
    factor::ComplexF64
    operators::Array{CWOperator}
end

copy(operator::CWLazyProduct) = CWLazyProduct(operator.factor,operator.operators)


function *(operator1::CWOperator,operator2::CWOperator)
    cpy1 = copy(operator1)
    cpy2 = copy(operator2)
    cpy1.factor = 1
    cpy2.factor = 1
    CWLazyProduct(operator1.factor*operator1.factor,[cpy1 cpy2])
end

function mul!(operator::CWOperator,state::cwstate)
    operator*state
end

function mul!(operator::CWLazyProduct,state::cwstate)
    operator.operators[end].factor = operator.operators[1].factor*operator.factor
    mul!(operator.operators[end],state::cwstate)
    operator.operators[end].factor = operator.operators[1].factor/operator.factor
    set_equal_diff!(state)
    for O in operator.operators[end-1:-1:2]
        mul!(O,state::cwstate)
        set_equal_diff!(state)
    end
    mul!(operator.operators[1],state::cwstate)       
end

function *(operator::CWLazyProduct,state::cwstate)
    save_state!(state)
    save_diff!(state)
    zero_diff!(state)
    operator.operators[end].factor = operator.operators[1].factor*operator.factor
    mul!(operator.operators[end],state::cwstate)
    operator.operators[end].factor = operator.operators[1].factor/operator.factor
    set_equal_diff!(state)
    for O in operator.operators[end-1:-1:2]
        mul!(O,state::cwstate)
        set_equal_diff!(state)
    end
    mul!(operator.operators[1],state::cwstate)
    load_diff!(state)
    load_state!(state)
end


function annihilation(b::cwbasis)
    return CWSingleOperator(1,annihilation!)
end
function creation(b::cwbasis)
    return CWSingleOperator(1,creation!)
end
function adw(b::cwbasis)
    return CWSingleOperator(1,adw!)
end
function wda(b::cwbasis)
    return CWSingleOperator(1,wda!)
end
function number(b::cwbasis)
    return CWSingleOperator(1,number!)
end

annihilation(state::cwstate) = annihilation(state.basis)
creation(state::cwstate) = creation(state.basis)
adw(state::cwstate) = adw(state.basis)
wda(state::cwstate) = wda(state.basis)
number(state::cwstate) = number(state.basis)

