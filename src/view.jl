"""
    TwoPhotonTimestepView{T}

Structure for viewing slice along same times in twophoton state.
"""
struct TwoPhotonTimestepView{T} <:AbstractVector{T}
    data::T
    timeindex::Int
    nsteps::Int
    offset::Int
end
function Base.getindex(x::TwoPhotonTimestepView,i::Int)
    if i<x.timeindex
        x.data[x.offset + twophoton_index(i,x.nsteps,x.timeindex)]
    #XXX Double Check!
    elseif i == x.timeindex
        sqrt(2)*x.data[x.offset + twophoton_index(x.timeindex,x.nsteps,x.timeindex)]
    else
        x.data[x.offset + twophoton_index(x.timeindex,x.nsteps,i)]
    end
end
Base.eachindex(x::TwoPhotonTimestepView) = 1:x.nsteps
function Base.setindex!(x::TwoPhotonTimestepView,input,i::Int)
    if i<x.timeindex
        x.data[x.offset + twophoton_index(i,x.nsteps,x.timeindex)] = input
    #XXX Double Check!
    elseif i == x.timeindex
        x.data[x.offset + twophoton_index(x.timeindex,x.nsteps,x.timeindex)]  = sqrt(2)*input
    else
        x.data[x.offset + twophoton_index(x.timeindex,x.nsteps,i)] = input
    end
end
Base.size(x::TwoPhotonTimestepView) = (x.nsteps,)



struct OnePhotonView{T} <: AbstractVector{T}
    data::AbstractArray{T}
    nsteps::Int
    offset::Int
end
Base.size(x::OnePhotonView) = (x.nsteps,)
function Base.getindex(x::OnePhotonView,i::Int)
    x.data[x.offset + i]
end
function Base.setindex!(x::OnePhotonView,input,i::Int)
    x.data[x.offset + i] = input
end

"""
    OnePhotonView(ψ::Ket)
    OnePhotonView(ψ::Ket,type)
    OnePhotonView(ψ::Ket,index)
    OnePhotonView(ψ::Ket,index,typw)


Return a view of the onephoton mode ``ξ(t)`` given a state defined on a [`WaveguideBasis`](@ref) or [`InputOutputWaveguideBasis`](@ref).
If the state is a [`InputOutputWaveguideBasis`](@ref) the `type` parameter can be used to choose between the input or output mode with `type = :input` or `type = :output`   
If no index is provided the ground state is returned. Index should follow same form outlined in [`view_waveguide`](@ref).
"""
function OnePhotonView(ψ::Ket;type=:input)
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    OnePhotonView(ψ,index,type=type)
end
function OnePhotonView(ψ::Ket,index;type=:input)
    viewed_data = view_waveguide(ψ::Ket,index)
    nsteps = get_nsteps(ψ.basis)
    if type==:input
        OnePhotonView(viewed_data,nsteps,1)
    elseif type==:output
        OnePhotonView(viewed_data,nsteps,1+nsteps)
    else
        error("Input type not recognized in OnePhotonView")
    end
end

"""
    TwophotonView{T} <: AbstractMatrix{T}

Structure for viewing twophoton state as symmetric matrix (only upper triangluar part is stored in memory).
"""
struct TwoPhotonView{T} <: AbstractMatrix{T}
    data::AbstractArray{T}
    nsteps::Int
    offset::Int
end
Base.size(x::TwoPhotonView) = (x.nsteps,x.nsteps)
function Base.getindex(x::TwoPhotonView,i::Int,j::Int)
    if i<j
        1/sqrt(2)*x.data[x.offset + twophoton_index(i,x.nsteps,j)]
    #XXX Double Check!
    elseif i == j
        x.data[x.offset + twophoton_index(i,x.nsteps,j)]
    else
        1/sqrt(2)*x.data[x.offset + twophoton_index(j,x.nsteps,i)]
    end
end
function Base.setindex!(x::TwoPhotonView,input,i::Int,j::Int)
    if i<j
        x.data[x.offset + twophoton_index(i,x.nsteps,j)] = sqrt(2)*input
    #XXX Double Check!
    elseif i == j
        x.data[x.offset + twophoton_index(i,x.nsteps,j)] = input
    else
        println("i>j")
        x.data[x.offset + twophoton_index(j,x.nsteps,i)] = sqrt(2)*input
    end
end
function twophoton_index(j,nsteps,timeindex)
    (j-1)*nsteps-((j-2)*(j-1))÷2+timeindex-j+1
end


function TwoPhotonView(ψ::Ket;type=:singlewaveguide)
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ::Ket,index;type=type)
end
function TwoPhotonView(ψ::Ket,index;type=:singlewaveguide)
    if type==:singlewaveguide
        viewed_data = view_waveguide(ψ::Ket,index)
        nsteps = get_nsteps(ψ.basis)
        TwoPhotonView(viewed_data,nsteps,nsteps+1)
    elseif type==:input
        viewed_data = view_waveguide(ψ::Ket,index)
        nsteps = get_nsteps(ψ.basis)
        TwoPhotonView(viewed_data,nsteps,2*nsteps+1)
    elseif type==:output
        viewed_data = view_waveguide(ψ::Ket,index)
        nsteps = get_nsteps(ψ.basis)
        TwoPhotonView(viewed_data,nsteps,2*nsteps+1+(nsteps*(nsteps+1))÷2)
    elseif type==:inputoutput
        viewed_data = view_waveguide(ψ::Ket,index)
        nsteps = get_nsteps(ψ.basis)
        InputOutputView(viewed_data,nsteps,2*nsteps+1+(nsteps*(nsteps+1)))
    else
        error("Input type not recognized in TwoPhotonView")
    end
end



"""
    InputOutputTimestepView{T} <:AbstractVector{T}

Structure for viewing slice along same times in one photon input output state.
"""
struct InputOutputTimestepView{T,I} <:AbstractVector{T}
    data::T
    timeindex::Int
    nsteps::Int
    offset::Int        
end
function InputOutputTimestepView(data::T,timeindex,nsteps,offset,type) where T
    InputOutputTimestepView{T,type}(data,timeindex,nsteps,offset)
end
function Base.getindex(x::InputOutputTimestepView{T,:input},i::Int) where T
    x.data[x.offset + (i-1)*x.nsteps + x.timeindex]
end
Base.eachindex(x::InputOutputTimestepView) = 1:x.nsteps
function Base.setindex!(x::InputOutputTimestepView{T,:input},input,i::Int) where T
    x.data[x.offset + (i-1)*x.nsteps + x.timeindex] = input
end
function Base.getindex(x::InputOutputTimestepView{T,:output},i::Int) where T
    x.data[x.offset + (x.timeindex-1)*x.nsteps + i]
end
function Base.setindex!(x::InputOutputTimestepView{T,:output},input,i::Int) where T
    x.data[x.offset + (x.timeindex-1)*x.nsteps + i] = input
end
Base.size(x::InputOutputTimestepView) = (x.nsteps,)


"""
InputOutputView{T} <: AbstractMatrix{T}

Structure for viewing state with one photon in input and output waveguide. 
"""
struct InputOutputView{T} <: AbstractMatrix{T}
    data::AbstractArray{T}
    nsteps::Int
    offset::Int
end
Base.size(x::InputOutputView) = (x.nsteps,x.nsteps)
function Base.getindex(x::InputOutputView,i::Int,j::Int)
    x.data[x.offset + (i-1)*x.nsteps + j]
end
function Base.setindex!(x::InputOutputView,input,i::Int,j::Int)
    x.data[x.offset + (i-1)*x.nsteps + j] = input
end
