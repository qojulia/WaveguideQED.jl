using QuantumOptics
using CavityWaveguide

"""
    InputOutputWaveguideBasis(N,times)

Basis for two waveguides containing in total `N` excitations. The two waveguides are refered to as input and output, but they are identical.
Currently restricted to either `N=1` or `N=2`. Times is timeinterval over which the photon state should be binned.
"""
mutable struct InputOutputWaveguideBasis{P} <: QuantumOptics.Basis
    shape::Vector{Int}
    N::Int
    offset::Int
    nsteps::Int
    function InputOutputWaveguideBasis(N,times)
        dim = 0
        if N==1
            dim = 1+2*length(times)
        elseif N==2
            dim = 1+2*length(times)+ length(times)^2 + 2*(length(times)^2 - (length(times)^2-length(times))÷2)
        else
            error("WCurrently N larger than 2 is not supported in InputOutputWaveguideBasis")
        end
        new{N}([dim+1], dim, 0,length(times))
    end
end

Base.:(==)(b1::InputOutputWaveguideBasis,b2::InputOutputWaveguideBasis) = (b1.N==b2.N && b1.offset==b2.offset && b1.nsteps==b2.nsteps)

function get_nsteps(basis::InputOutputWaveguideBasis)
    basis.nsteps
end

"""
    zerophoton(bw::InputOutputWaveguideBasis)

Create a waveguide vacuum state ``|\\emptyset,\\emptyset \\rangle`` of the two waveguides.
"""
function zerophoton(b::InputOutputWaveguideBasis)
    state = Ket(b)
    state.data[1] = 1
    return state
end


"""
    onephoton(b::InputOutputWaveguideBasis,type,ξ::Function,times,args...,norm=True)
    onephoton(b::InputOutputWaveguideBasis,type,ξ;norm=true)

Create a onephoton wavepacket of the form ``W_{i/o}^\\dagger(\\xi) |\\emptyset,\\emptyset \\rangle = \\int_{t_0}^{t_{end}} dt  \\xi(t) w_{i/o}^\\dagger(t) |\\emptyset,\\emptyset \\rangle`` where the subscript i/o denotes the input or output waveguide and determined by the `type` parameter.

### Input
* `b` is the basis of the input output waveguides.    
*  `type` determines which waveguide the onephoton state is created in. Can be either `:input` or `:outputput`.
* ξ can either be a function evaluated as `ξ.(times,args...)` or a vector of length: `b.nsteps`.
* If `norm==true` the state is normalized through `normalize!`.
"""
function onephoton(basis::InputOutputWaveguideBasis,type,ξ::Function,times,args...;norm=true)
    onephoton(basis,type,ξ.(times,args...),norm=norm)
end
function onephoton(basis::InputOutputWaveguideBasis,type,ξvec;norm=true)
    state = Ket(basis)
    view = OnePhotonView(state,type=type)
    view .= ξvec
    if norm
        normalize!(state)
    end
    return state
end

"""
    twophoton(b::InputOutputWaveguideBasis,type,ξ::Function,times,args...,norm=True)
    twophoton(b::InputOutputWaveguideBasis,type,ξ;norm=true)

Create a twophoton wavepacket of the form ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{i/o}^\\dagger(t)w_{i/o}^\\dagger(t') |\\emptyset,\\emptyset \\rangle``
The subscript i/o denotes the input or output waveguide and is determined by the `type` parameter.

### Input
* `b` is the basis of the input output waveguides.    
*  `type` determines which waveguide the twophoton state is created in. Can be either `:input` or `:outputput` or `:inputoutput` for an entangled onephoton state in the input waveguide and onephoton state in the outputwaveguide.
* ξ can either be a function evaluated as `ξ.(times,args...)` or a vector of length: `b.nsteps`.
* If `norm==true` the state is normalized through `normalize!`.
"""
function twophoton(b::InputOutputWaveguideBasis,type,ξ::Function,times,args...;norm=true)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,type=type)
    if type==:inputoutput
        for i in 1:nsteps
            for j in 1:nsteps
                viewed_data[i,j] = ξ(times[i],times[j],args...)
            end
        end    
    else
        for i in 1:nsteps
            for j in i:nsteps
                viewed_data[i,j] = ξ(times[i],times[j],args...)
            end
        end
    end
    if norm
        normalize!(state)
    end
    return state
end
function twophoton(b::InputOutputWaveguideBasis,type,ξ::Matrix;norm=true)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,type=type)
    if type==:inputoutput
        for i in 1:nsteps
            for j in 1:nsteps
                viewed_data[i,j] = ξ[i,j]
            end
        end    
    else
        for i in 1:nsteps
            for j in i:nsteps
                viewed_data[i,j] = ξ[i,j]
            end
        end
    end
    if norm
        normalize!(state)
    end
    return state
end

get_waveguide_location(basis::InputOutputWaveguideBasis) = 1

mutable struct InputWaveguideCreate{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end
mutable struct InputWaveguideDestroy{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end
mutable struct OutputWaveguideCreate{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end
mutable struct OutputWaveguideDestroy{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end

"""
    dagger(op::WaveguideCreate)
    dagger(op::WaveguideDestroy)

Dagger opration on Waveguide operator. 

"""
function dagger(op::InputWaveguideCreate)
    @assert op.basis_l == op.basis_r
    out = inputdestroy(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out 
end
function dagger(op::InputWaveguideDestroy)
    @assert op.basis_l == op.basis_r
    out = inputcreate(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out
end
function dagger(op::OutputWaveguideCreate)
    @assert op.basis_l == op.basis_r
    out = outputdestroy(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out 
end
function dagger(op::OutputWaveguideDestroy)
    @assert op.basis_l == op.basis_r
    out = outputcreate(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out
end


function Base.:copy(x::InputWaveguideDestroy{B,B,N}) where {B,N}
    InputWaveguideDestroy{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end
function Base.:copy(x::OutputWaveguideDestroy{B,B,N}) where {B,N}
    OutputWaveguideDestroy{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end
function Base.:copy(x::InputWaveguideCreate{B,B,N}) where {B,N}
    InputWaveguideCreate{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end
function Base.:copy(x::OutputWaveguideCreate{B,B,N}) where {B,N}
    OutputWaveguideCreate{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end


#Destroy 1 input waveguide photon
function waveguide_mul!(result,a::InputWaveguideDestroy{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += a.factor*alpha*b[1+timeindex]
    return result
end
#Destroy 1 input waveguide photon
function waveguide_mul!(result,a::OutputWaveguideDestroy{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += a.factor*alpha*b[1+nsteps+timeindex]
    return result
end
#Create 1 input waveguide photon
function waveguide_mul!(result,a::InputWaveguideCreate{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1+timeindex] += a.factor*alpha*b[1]
    return result
end
#Create 1 input waveguide photon
function waveguide_mul!(result,a::OutputWaveguideCreate{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1+nsteps+timeindex] += a.factor*alpha*b[1]
    return result
end

#Destroy 2 input waveguide photon
function waveguide_mul!(result,a::InputWaveguideDestroy{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += a.factor*alpha*b[1+timeindex]
    twophotonview_input = TwoPhotonTimestepView(b,timeindex,nsteps,1+2*nsteps)
    twophotonview_input_output = InputOutputTimestepView(b,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:output)
    axpy!(alpha*a.factor,twophotonview_input_output,view(result,2+nsteps:1:2*nsteps+1))
    axpy!(alpha*a.factor,twophotonview_input,view(result,2:1:nsteps+1))
    return result
end

#Destroy 2 output waveguide photon
function waveguide_mul!(result,a::OutputWaveguideDestroy{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += a.factor*alpha*b[1+nsteps+timeindex]
    twophotonview_output = TwoPhotonTimestepView(b,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1))÷2)
    twophotonview_input_output = InputOutputTimestepView(b,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:input)
    axpy!(alpha*a.factor,twophotonview_output,view(result,2+nsteps:1:2*nsteps+1))
    axpy!(alpha*a.factor,twophotonview_input_output,view(result,2:1:nsteps+1))
    return result
end

#Create 2 waveguide photon
function waveguide_mul!(result,a::InputWaveguideCreate{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1+timeindex] += a.factor*alpha*b[1]
    twophotonview_input = TwoPhotonTimestepView(result,timeindex,nsteps,1+2*nsteps)
    twophotonview_input_output = InputOutputTimestepView(result,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:output)
    axpy!(alpha*a.factor,view(b,2+nsteps:1:2*nsteps+1),twophotonview_input_output)
    axpy!(alpha*a.factor,view(b,2:1:nsteps+1),twophotonview_input)
    return result
end

#Destroy 2 output waveguide photon
function waveguide_mul!(result,a::OutputWaveguideCreate{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds  result[1+nsteps+timeindex] +=  a.factor*alpha*b[1]
    twophotonview_output = TwoPhotonTimestepView(result,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1))÷2)
    twophotonview_input_output = InputOutputTimestepView(result,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:input)
    axpy!(alpha*a.factor,view(b,2+nsteps:1:2*nsteps+1),twophotonview_output)
    axpy!(alpha*a.factor,view(b,2:1:nsteps+1),twophotonview_input_output)
    return result
end

"""
    inputdestroy(basis::WaveguideBasis{N})

Annihilation operator for input waveguide in [`InputOutputWaveguideBasis`](@ref) for either one or two photons. 

"""
function inputdestroy(basis::InputOutputWaveguideBasis{N}) where N
    B = typeof(basis)
    return InputWaveguideDestroy{B,B,N}(basis,basis,1,1)
end

"""
    inputcreate(basis::WaveguideBasis{N}) where N

Creation operator for input waveguide in [`InputOutputWaveguideBasis`](@ref) for either one or two photons. 

"""
function inputcreate(basis::InputOutputWaveguideBasis{N}) where N
    B = typeof(basis)
    return InputWaveguideCreate{B,B,N}(basis,basis,1,1)
end


"""
    outputdestroy(basis::WaveguideBasis{N})

Annihilation operator for output waveguide in [`InputOutputWaveguideBasis`](@ref) for either one or two photons. 

"""
function outputdestroy(basis::InputOutputWaveguideBasis{N}) where N
    B = typeof(basis)
    return OutputWaveguideDestroy{B,B,N}(basis,basis,1,1)
end

"""
    Outputcreate(basis::WaveguideBasis{N}) where N

Creation operator for output waveguide in [`InputOutputWaveguideBasis`](@ref) for either one or two photons. 
"""
function outputcreate(basis::InputOutputWaveguideBasis{N}) where N
    B = typeof(basis)
    return OutputWaveguideCreate{B,B,N}(basis,basis,1,1)
end




"""
    inputabsorption(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    inputabsorption(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T

Create [`CavityWaveguideAbsorption`](@ref) that applies `create(b::FockBasis)` on `FockBasis` and inputdestroy(b::InputOutputWaveguideBasis{T}) on [`InputOutputWaveguideBasis{T}`](@ref).  
"""
function inputabsorption(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),inputdestroy(b1),[1,2])
end
function inputabsorption(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),inputdestroy(b2),[2,1])
end

"""
    inputemission(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    inputemission(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T

Create [`CavityWaveguideEmission`](@ref) that applies `destroy(b::FockBasis)` on `FockBasis` and inputcreate(b::InputOutputWaveguideBasis{T}) on [`InputOutputWaveguideBasis{T}`](@ref).  
"""
function inputemission(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),inputcreate(b1),[1,2])
end
function inputemission(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),inputcreate(b2),[2,1])
end

"""
    outputabsorption(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    outputabsorption(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T

Create [`CavityWaveguideAbsorption`](@ref) that applies `create(b::FockBasis)` on `FockBasis` and outputdestroy(b::InputOutputWaveguideBasis{T}) on [`InputOutputWaveguideBasis{T}`](@ref).  
"""
function outputabsorption(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),outputdestroy(b1),[1,2])
end
function outputabsorption(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),outputdestroy(b2),[2,1])
end

"""
    outputemission(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    outputemission(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T

Create [`CavityWaveguideEmission`](@ref) that applies `destroy(b::FockBasis)` on `FockBasis` and outputcreate(b::InputOutputWaveguideBasis{T}) on [`InputOutputWaveguideBasis{T}`](@ref).  
"""
function outputemission(b1::InputOutputWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),outputcreate(b1),[1,2])
end
function outputemission(b1::FockBasis,b2::InputOutputWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),outputcreate(b2),[2,1])
end

Base.:*(x::WaveguideOperator{B1,B2},y::WaveguideOperator{B1,B2}) where {B1,B2} = LazyProduct(x,y)
