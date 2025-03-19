"""
    WaveguideDispersoin{B1,B2,N,idx} <: WaveguideOperator{B1,B2}

Operator structure for dispatching dispersion operation ``\\sum_n w_n^\\dagger w_{n-1} + w_{n-1}^\\dagger w_n ``.
Np is used to dispatch one or two photon routine and idx denotes the index of the waveguide the operator is acting on. 
"""
mutable struct WaveguideDispersion{B1,B2,N,idx} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
end


function Base.:copy(x::WaveguideDispersion{B,B,Np,idx}) where {B,Np,idx}
    WaveguideDispersion{B,B,Np,idx}(x.basis_l,x.basis_r,x.factor)
end

@inline function set_time!(o::WaveguideDispersion, t::Number)
    return o
end
function dispersion(basis::WaveguideBasis{Np,1}) where {Np}
    B = typeof(basis)
    return WaveguideDispersion{B,B,Np,1}(basis,basis,1)
end
function dispersion(basis::WaveguideBasis{Np,Nw},i::Int) where {Np,Nw}
    @assert i <= Nw
    B = typeof(basis)
    return WaveguideDispersion{B,B,Np,i}(basis,basis,1)
end


function waveguide_mul!(result,a::WaveguideDispersion{B,B,1,idx},b,alpha,beta) where {B,idx}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    nsteps = a.basis_l.nsteps
    
    #@inbounds result[2+(idx-1)*nsteps:(idx)*nsteps+1] .+= 1/0.01*b[2+(idx-1)*nsteps:(idx)*nsteps+1]
    
    @inbounds result[2+(idx-1)*nsteps:(idx)*nsteps] .+= alpha*a.factor*b[3+(idx-1)*nsteps:1+(idx)*nsteps]
    @inbounds result[3+(idx-1)*nsteps:1+(idx)*nsteps] .+= alpha*a.factor*b[2+(idx-1)*nsteps:(idx)*nsteps]
    #@inbounds result[2+(idx-1)*nsteps:1+(idx)*nsteps] .+= alpha*a.factor*b[1+(idx-1)*nsteps:(idx)*nsteps]
    return
end