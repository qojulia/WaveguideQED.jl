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
function Base.setindex!(x::TwoPhotonTimestepView,Left,i::Int)
    if i<x.timeindex
        x.data[x.offset + twophoton_index(i,x.nsteps,x.timeindex)] = Left
    #XXX Double Check!
    elseif i == x.timeindex
        x.data[x.offset + twophoton_index(x.timeindex,x.nsteps,x.timeindex)]  = sqrt(2)*Left
    else
        x.data[x.offset + twophoton_index(x.timeindex,x.nsteps,i)] = Left
    end
end
Base.size(x::TwoPhotonTimestepView) = (x.nsteps,)
function axpy!(α, x::AbstractArray, y::TwoPhotonTimestepView)
    n = length(x)
    if n != length(y)
        throw(DimensionMismatch("x has length $n, but y has length $(length(y))"))
    end
    for (IY, IX) in zip(eachindex(y), eachindex(x))
        if IY == y.timeindex
            @inbounds y[IY] = y[IY]/2+ x[IX]*α
        else
            @inbounds y[IY] += x[IX]*α
        end
    end
    y
end

struct OnePhotonView{T} <: AbstractVector{T}
    data::AbstractArray{T}
    nsteps::Int
    offset::Int
end
Base.size(x::OnePhotonView) = (x.nsteps,)
function Base.getindex(x::OnePhotonView,i::Int)
    x.data[x.offset + i]
end
function Base.setindex!(x::OnePhotonView,Left,i::Int)
    x.data[x.offset + i] = Left
end

"""
    OnePhotonView(ψ::Ket)
    OnePhotonView(ψ::Ket,type)
    OnePhotonView(ψ::Ket,index)
    OnePhotonView(ψ::Ket,index,typw)


Return a view of the onephoton mode ``ξ(t)`` given a state defined on a [`WaveguideBasis`](@ref) or [`LeftRightWaveguideBasis`](@ref).
If the state is a [`LeftRightWaveguideBasis`](@ref) the `type` parameter can be used to choose between the Left or Right mode with `type = :Left` or `type = :Right`   
If no index is provided the ground state is returned. Index should follow same form outlined in [`view_waveguide`](@ref).
"""
function OnePhotonView(ψ::Ket;type=:Left)
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    OnePhotonView(ψ,index,type=type)
end
function OnePhotonView(ψ::Ket,index;type=:Left)
    viewed_data = view_waveguide(ψ::Ket,index)
    nsteps = get_nsteps(ψ.basis)
    if type==:Left
        OnePhotonView(viewed_data,nsteps,1)
    elseif type==:Right
        OnePhotonView(viewed_data,nsteps,1+nsteps)
    else
        error("Left type not recognized in OnePhotonView")
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
function Base.setindex!(x::TwoPhotonView,Left,i::Int,j::Int)
    if i<j
        x.data[x.offset + twophoton_index(i,x.nsteps,j)] = sqrt(2)*Left
    #XXX Double Check!
    elseif i == j
        x.data[x.offset + twophoton_index(i,x.nsteps,j)] = Left
    else
        println("i>j")
        x.data[x.offset + twophoton_index(j,x.nsteps,i)] = sqrt(2)*Left
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
    elseif type==:Left
        viewed_data = view_waveguide(ψ::Ket,index)
        nsteps = get_nsteps(ψ.basis)
        TwoPhotonView(viewed_data,nsteps,2*nsteps+1)
    elseif type==:Right
        viewed_data = view_waveguide(ψ::Ket,index)
        nsteps = get_nsteps(ψ.basis)
        TwoPhotonView(viewed_data,nsteps,2*nsteps+1+(nsteps*(nsteps+1))÷2)
    elseif type==:LeftRight
        viewed_data = view_waveguide(ψ::Ket,index)
        nsteps = get_nsteps(ψ.basis)
        LeftRightView(viewed_data,nsteps,2*nsteps+1+(nsteps*(nsteps+1)))
    else
        error("Left type not recognized in TwoPhotonView")
    end
end



"""
    LeftRightTimestepView{T} <:AbstractVector{T}

Structure for viewing slice along same times in one photon Left Right state.
"""
struct LeftRightTimestepView{T,I} <:AbstractVector{T}
    data::T
    timeindex::Int
    nsteps::Int
    offset::Int        
end
function LeftRightTimestepView(data::T,timeindex,nsteps,offset,type) where T
    LeftRightTimestepView{T,type}(data,timeindex,nsteps,offset)
end
function Base.getindex(x::LeftRightTimestepView{T,:Left},i::Int) where T
    x.data[x.offset + (i-1)*x.nsteps + x.timeindex]
end
Base.eachindex(x::LeftRightTimestepView) = 1:x.nsteps
function Base.setindex!(x::LeftRightTimestepView{T,:Left},Left,i::Int) where T
    x.data[x.offset + (i-1)*x.nsteps + x.timeindex] = Left
end
function Base.getindex(x::LeftRightTimestepView{T,:Right},i::Int) where T
    x.data[x.offset + (x.timeindex-1)*x.nsteps + i]
end
function Base.setindex!(x::LeftRightTimestepView{T,:Right},Left,i::Int) where T
    x.data[x.offset + (x.timeindex-1)*x.nsteps + i] = Left
end
Base.size(x::LeftRightTimestepView) = (x.nsteps,)


"""
LeftRightView{T} <: AbstractMatrix{T}

Structure for viewing state with one photon in Left and Right waveguide. 
"""
struct LeftRightView{T} <: AbstractMatrix{T}
    data::AbstractArray{T}
    nsteps::Int
    offset::Int
end
Base.size(x::LeftRightView) = (x.nsteps,x.nsteps)
function Base.getindex(x::LeftRightView,i::Int,j::Int)
    x.data[x.offset + (i-1)*x.nsteps + j]
end
function Base.setindex!(x::LeftRightView,Left,i::Int,j::Int)
    x.data[x.offset + (i-1)*x.nsteps + j] = Left
end
