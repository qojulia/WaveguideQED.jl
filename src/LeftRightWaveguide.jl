using QuantumOptics
using CavityWaveguide

"""
    LeftRightWaveguideBasis(N,times)

Basis for two waveguides containing in total `N` excitations. The two waveguides are refered to as Left and Right, but they are identical.
Currently restricted to either `N=1` or `N=2`. Times is timeinterval over which the photon state should be binned.
"""
mutable struct LeftRightWaveguideBasis{P} <: QuantumOptics.Basis
    shape::Vector{Int}
    N::Int
    offset::Int
    nsteps::Int
    function LeftRightWaveguideBasis(N,times)
        dim = 0
        if N==1
            dim = 1+2*length(times)
        elseif N==2
            dim = 1+2*length(times)+ length(times)^2 + 2*(length(times)^2 - (length(times)^2-length(times))÷2)
        else
            error("WCurrently N larger than 2 is not supported in LeftRightWaveguideBasis")
        end
        new{N}([dim+1], dim, 0,length(times))
    end
end

Base.:(==)(b1::LeftRightWaveguideBasis,b2::LeftRightWaveguideBasis) = (b1.N==b2.N && b1.offset==b2.offset && b1.nsteps==b2.nsteps)

function get_nsteps(basis::LeftRightWaveguideBasis)
    basis.nsteps
end

"""
    zerophoton(bw::LeftRightWaveguideBasis)

Create a waveguide vacuum state ``|\\emptyset,\\emptyset \\rangle`` of the two waveguides.
"""
function zerophoton(b::LeftRightWaveguideBasis)
    state = Ket(b)
    state.data[1] = 1
    return state
end


"""
    onephoton(b::LeftRightWaveguideBasis,type,ξ::Function,times,args...,norm=True)
    onephoton(b::LeftRightWaveguideBasis,type,ξ;norm=true)

Create a onephoton wavepacket of the form ``W_{i/o}^\\dagger(\\xi) |\\emptyset,\\emptyset \\rangle = \\int_{t_0}^{t_{end}} dt  \\xi(t) w_{i/o}^\\dagger(t) |\\emptyset,\\emptyset \\rangle`` where the subscript i/o denotes the Left or Right waveguide and determined by the `type` parameter.

### Input
* `b` is the basis of the Left Right waveguides.    
*  `type` determines which waveguide the onephoton state is created in. Can be either `:Left` or `:Right`.
* ξ can either be a function evaluated as `ξ.(times,args...)` or a vector of length: `b.nsteps`.
* If `norm==true` the state is normalized through `normalize!`.
"""
function onephoton(basis::LeftRightWaveguideBasis,type,ξ::Function,times,args...;norm=true)
    onephoton(basis,type,ξ.(times,args...),norm=norm)
end
function onephoton(basis::LeftRightWaveguideBasis,type,ξvec;norm=true)
    state = Ket(basis)
    view = OnePhotonView(state,type=type)
    view .= ξvec
    if norm
        normalize!(state)
    end
    return state
end

"""
    twophoton(b::LeftRightWaveguideBasis,type,ξ::Function,times,args...,norm=True)
    twophoton(b::LeftRightWaveguideBasis,type,ξ;norm=true)

Create a twophoton wavepacket of the form ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{i/o}^\\dagger(t)w_{i/o}^\\dagger(t') |\\emptyset,\\emptyset \\rangle``
The subscript i/o denotes the Left or Right waveguide and is determined by the `type` parameter.

### Input
* `b` is the basis of the Left Right waveguides.    
*  `type` determines which waveguide the twophoton state is created in. Can be either `:Left` or `:Right` or `:LeftRight` for an entangled onephoton state in the Left waveguide and onephoton state in the Rightwaveguide.
* ξ can either be a function evaluated as `ξ.(times,args...)` or a vector of length: `b.nsteps`.
* If `norm==true` the state is normalized through `normalize!`.
"""
function twophoton(b::LeftRightWaveguideBasis,type,ξ::Function,times,args...;norm=true)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,type=type)
    if type==:LeftRight
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
function twophoton(b::LeftRightWaveguideBasis,type,ξ::Matrix;norm=true)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,type=type)
    if type==:LeftRight
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

get_waveguide_location(basis::LeftRightWaveguideBasis) = 1

mutable struct LeftWaveguideCreate{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end
mutable struct LeftWaveguideDestroy{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end
mutable struct RightWaveguideCreate{B1,B2,N} <:WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    timeindex::Int
end
mutable struct RightWaveguideDestroy{B1,B2,N} <:WaveguideOperator{B1,B2}
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
function dagger(op::LeftWaveguideCreate)
    @assert op.basis_l == op.basis_r
    out = Leftdestroy(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out 
end
function dagger(op::LeftWaveguideDestroy)
    @assert op.basis_l == op.basis_r
    out = Leftcreate(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out
end
function dagger(op::RightWaveguideCreate)
    @assert op.basis_l == op.basis_r
    out = Rightdestroy(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out 
end
function dagger(op::RightWaveguideDestroy)
    @assert op.basis_l == op.basis_r
    out = Rightcreate(op.basis_l)
    out.factor = op.factor
    out.timeindex = op.timeindex
    out
end


function Base.:copy(x::LeftWaveguideDestroy{B,B,N}) where {B,N}
    LeftWaveguideDestroy{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end
function Base.:copy(x::RightWaveguideDestroy{B,B,N}) where {B,N}
    RightWaveguideDestroy{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end
function Base.:copy(x::LeftWaveguideCreate{B,B,N}) where {B,N}
    LeftWaveguideCreate{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end
function Base.:copy(x::RightWaveguideCreate{B,B,N}) where {B,N}
    RightWaveguideCreate{B,B,N}(x.basis_l,x.basis_r,x.factor,x.timeindex)
end


#Destroy 1 Left waveguide photon
function waveguide_mul!(result,a::LeftWaveguideDestroy{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    @inbounds result[1] += a.factor*alpha*b[1+timeindex]
    @assert sum(isnan.(result)) == 0
    return result
end
#Destroy 1 Left waveguide photon
function waveguide_mul!(result,a::RightWaveguideDestroy{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += a.factor*alpha*b[1+nsteps+timeindex]
    @assert sum(isnan.(result)) == 0
    return result
end
#Create 1 Left waveguide photon
function waveguide_mul!(result,a::LeftWaveguideCreate{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    @inbounds result[1+timeindex] += a.factor*alpha*b[1]
    @assert sum(isnan.(result)) == 0
    return result
end
#Create 1 Left waveguide photon
function waveguide_mul!(result,a::RightWaveguideCreate{B,B,1},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1+nsteps+timeindex] += a.factor*alpha*b[1]
    return result
end

#Destroy 2 Left waveguide photon
function waveguide_mul!(result,a::LeftWaveguideDestroy{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += a.factor*alpha*b[1+timeindex]
    twophotonview_Left = TwoPhotonTimestepView(b,timeindex,nsteps,1+2*nsteps)
    twophotonview_Left_Right = LeftRightTimestepView(b,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:Right)
    axpy!(alpha*a.factor,twophotonview_Left_Right,view(result,2+nsteps:1:2*nsteps+1))
    axpy!(alpha*a.factor,twophotonview_Left,view(result,2:1:nsteps+1))
    #@assert sum(isnan.(result)) == 0
    return result
end

#Destroy 2 Right waveguide photon
function waveguide_mul!(result,a::RightWaveguideDestroy{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1] += a.factor*alpha*b[1+nsteps+timeindex]
    twophotonview_Right = TwoPhotonTimestepView(b,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1))÷2)
    twophotonview_Left_Right = LeftRightTimestepView(b,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:Left)
    axpy!(alpha*a.factor,twophotonview_Right,view(result,2+nsteps:1:2*nsteps+1))
    axpy!(alpha*a.factor,twophotonview_Left_Right,view(result,2:1:nsteps+1))
    #@assert sum(isnan.(result)) == 0
    
    return result
end

#Create 2 waveguide photon
function waveguide_mul!(result,a::LeftWaveguideCreate{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds result[1+timeindex] += a.factor*alpha*b[1]
    twophotonview_Left = TwoPhotonTimestepView(result,timeindex,nsteps,1+2*nsteps)
    twophotonview_Left_Right = LeftRightTimestepView(result,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:Right)
    axpy!(alpha*a.factor,view(b,2+nsteps:1:2*nsteps+1),twophotonview_Left_Right)
    axpy!(alpha*a.factor,view(b,2:1:nsteps+1),twophotonview_Left)
    #@assert sum(isnan.(result)) == 0
    
    return result
end

#Destroy 2 Right waveguide photon
function waveguide_mul!(result,a::RightWaveguideCreate{B,B,2},b,alpha,beta) where {B}
    if !isone(beta)
        rmul!(result,beta)
    end
    timeindex = a.timeindex
    nsteps = a.basis_l.nsteps
    @inbounds  result[1+nsteps+timeindex] +=  a.factor*alpha*b[1]
    twophotonview_Right = TwoPhotonTimestepView(result,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1))÷2)
    twophotonview_Left_Right = LeftRightTimestepView(result,timeindex,nsteps,1+2*nsteps+(nsteps*(nsteps+1)),:Left)
    axpy!(alpha*a.factor,view(b,2+nsteps:1:2*nsteps+1),twophotonview_Right)
    axpy!(alpha*a.factor,view(b,2:1:nsteps+1),twophotonview_Left_Right)
    #@assert sum(isnan.(result)) == 0
    
    return result
end

"""
    Leftdestroy(basis::WaveguideBasis{N})

Annihilation operator for Left waveguide in [`LeftRightWaveguideBasis`](@ref) for either one or two photons. 

"""
function Leftdestroy(basis::LeftRightWaveguideBasis{N}) where N
    B = typeof(basis)
    return LeftWaveguideDestroy{B,B,N}(basis,basis,1,1)
end

"""
    Leftcreate(basis::WaveguideBasis{N}) where N

Creation operator for Left waveguide in [`LeftRightWaveguideBasis`](@ref) for either one or two photons. 

"""
function Leftcreate(basis::LeftRightWaveguideBasis{N}) where N
    B = typeof(basis)
    return LeftWaveguideCreate{B,B,N}(basis,basis,1,1)
end


"""
    Rightdestroy(basis::WaveguideBasis{N})

Annihilation operator for Right waveguide in [`LeftRightWaveguideBasis`](@ref) for either one or two photons. 

"""
function Rightdestroy(basis::LeftRightWaveguideBasis{N}) where N
    B = typeof(basis)
    return RightWaveguideDestroy{B,B,N}(basis,basis,1,1)
end

"""
    Rightcreate(basis::WaveguideBasis{N}) where N

Creation operator for Right waveguide in [`LeftRightWaveguideBasis`](@ref) for either one or two photons. 
"""
function Rightcreate(basis::LeftRightWaveguideBasis{N}) where N
    B = typeof(basis)
    return RightWaveguideCreate{B,B,N}(basis,basis,1,1)
end




"""
    Leftabsorption(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    Leftabsorption(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T

Create [`CavityWaveguideAbsorption`](@ref) that applies `create(b::FockBasis)` on `FockBasis` and Leftdestroy(b::LeftRightWaveguideBasis{T}) on [`LeftRightWaveguideBasis{T}`](@ref).  
"""
function Leftabsorption(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),Leftdestroy(b1),[1,2])
end
function Leftabsorption(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),Leftdestroy(b2),[2,1])
end

"""
    Leftemission(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    Leftemission(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T

Create [`CavityWaveguideEmission`](@ref) that applies `destroy(b::FockBasis)` on `FockBasis` and Leftcreate(b::LeftRightWaveguideBasis{T}) on [`LeftRightWaveguideBasis{T}`](@ref).  
"""
function Leftemission(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),Leftcreate(b1),[1,2])
end
function Leftemission(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),Leftcreate(b2),[2,1])
end

"""
    Rightabsorption(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    Rightabsorption(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T

Create [`CavityWaveguideAbsorption`](@ref) that applies `create(b::FockBasis)` on `FockBasis` and Rightdestroy(b::LeftRightWaveguideBasis{T}) on [`LeftRightWaveguideBasis{T}`](@ref).  
"""
function Rightabsorption(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),Rightdestroy(b1),[1,2])
end
function Rightabsorption(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideAbsorption(btotal,btotal,complex(1.0),Rightdestroy(b2),[2,1])
end

"""
    Rightemission(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    Rightemission(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T

Create [`CavityWaveguideEmission`](@ref) that applies `destroy(b::FockBasis)` on `FockBasis` and Rightcreate(b::LeftRightWaveguideBasis{T}) on [`LeftRightWaveguideBasis{T}`](@ref).  
"""
function Rightemission(b1::LeftRightWaveguideBasis{T},b2::FockBasis) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),Rightcreate(b1),[1,2])
end
function Rightemission(b1::FockBasis,b2::LeftRightWaveguideBasis{T}) where T
    btotal = tensor(b1,b2)
    return CavityWaveguideEmission(btotal,btotal,complex(1.0),Rightcreate(b2),[2,1])
end

