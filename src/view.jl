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
function axpy!(α, x::TwoPhotonTimestepView, y::TwoPhotonTimestepView)
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
    OnePhotonView(ψ::T) where {T<:SingleWaveguideKet}
    OnePhotonView(ψ::T,index::I) where {T<:SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    OnePhotonView(ψ::T,WI::Int)  where {T<:MultipleWaveguideKet}
    OnePhotonView(ψ::T,index::I,WI::Int) where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}

Return a view of the onephoton wavefunction ``ξ(t)`` given a state containing a [`WaveguideBasis`](@ref).
If the [`WaveguideBasis`](@ref) contains more than one waveguide, a waveguide index ``WI`` is required to indicate which waveguide is viewed (1,2,3,... etc.)  
If the state contains more basises (e.g. a cavity) ``index`` used to indicate which state should be viewed. Index should follow same form outlined in [`view_waveguide`](@ref) and if not given the groundstate is assummed.
"""
function OnePhotonView(ψ::T) where {T<:SingleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    OnePhotonView(ψ,index)
end
function OnePhotonView(ψ::T,index::I) where {T<:SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    OnePhotonView(viewed_data,nsteps,1)
end
function OnePhotonView(ψ::T,WI::Int)  where {T<:MultipleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    OnePhotonView(ψ,index,WI)
end
function OnePhotonView(ψ::T,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    error("When considering multiple waveguides please provide also index for the waveguide")
end
function OnePhotonView(ψ::T,index::I,WI::Int) where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    @assert WI <= get_number_of_waveguides(ψ.basis)
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    OnePhotonView(viewed_data,nsteps,1+nsteps*(WI-1))
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

"""
    TwoPhotonView(ψ::T) where {T <: SingleWaveguideKet}
    TwoPhotonView(ψ::T,index::I) where {T <: SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    TwoPhotonView(ψ::T,WI::Int) where {T <: MultipleWaveguideKet}
    TwoPhotonView(ψ::T,WI1::Int,WI2::Int) where {T <: MultipleWaveguideKet}
    TwoPhotonView(ψ::T,WI::I) where {T <: MultipleWaveguideKet,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}}
    TwoPhotonView(ψ::T,index::I,WI::Int)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    
Return a view of the twophoton wavefunction ``ξ(t1,t2)`` given a state containing a [`WaveguideBasis`](@ref).
If the [`WaveguideBasis`](@ref) contains more than one waveguide, a waveguide index ``WI`` is required to indicate which waveguide is viewed.
`WI` follows the same syntax  outlined in [`twophoton`](@ref) for more information on how to view the state.
If the state contains more basises (e.g. a cavity) ``index`` used to indicate which state should be viewed. Index should follow same form outlined in [`view_waveguide`](@ref) and if not given the groundstate is assummed.
"""
function TwoPhotonView(ψ::T) where {T <: SingleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,index)
end
function TwoPhotonView(ψ::T,index::I) where {T <: SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    TwoPhotonView(viewed_data,nsteps,nsteps+1)
end
function TwoPhotonView(ψ::T,WI::Int) where {T <: MultipleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,index,WI)
end
function TwoPhotonView(ψ::T,WI1::Int,WI2::Int) where {T <: MultipleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,index,WI1,WI2)
end
function TwoPhotonView(ψ::T,WI::I) where {T <: MultipleWaveguideKet,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,index,WI)
end
function TwoPhotonView(ψ::T,index::I,WI::Int)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    Nw = get_number_of_waveguides(ψ.basis)
    @assert WI<=Nw
    TwoPhotonView(viewed_data,nsteps,Nw*nsteps+1+(WI-1)*(nsteps*(nsteps+1)÷2))
end
function TwoPhotonView(ψ::T,index::I,WI1::Int,WI2::Int)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    Nw = get_number_of_waveguides(ψ.basis)
    i,j = min(WI1,WI2),max(WI1,WI2)
    @assert i<=Nw && j<=Nw
    index = (i-1)*Nw + j - (i*(i+1))÷2
    TwoWaveguideView(viewed_data,nsteps,Nw*nsteps+1+(Nw)*(nsteps*(nsteps+1)÷2)+(index-1)*nsteps^2,i==WI1)
end
TwoPhotonView(ψ::T,index::I,WI::F)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}},F<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = TwoPhotonView(ψ,index,WI[1],WI[2]) 

"""
TwoWaveguideTimestepView{T} <:AbstractVector{T}

Structure for viewing slice along same times in one photon Left Right state.
"""
struct TwoWaveguideTimestepView{T,I} <:AbstractVector{T}
    data::T
    timeindex::Int
    nsteps::Int
    offset::Int        
end
function TwoWaveguideTimestepView(data::T,timeindex,nsteps,offset,indexing) where T
    TwoWaveguideTimestepView{T,indexing}(data,timeindex,nsteps,offset)
end
function Base.getindex(x::TwoWaveguideTimestepView{T,true},i::Int) where T
    x.data[x.offset + (i-1)*x.nsteps + x.timeindex]
end
function Base.setindex!(x::TwoWaveguideTimestepView{T,true},input,i::Int) where T
    x.data[x.offset + (i-1)*x.nsteps + x.timeindex] = input
end
function Base.getindex(x::TwoWaveguideTimestepView{T,false},i::Int) where T
    x.data[x.offset + (x.timeindex-1)*x.nsteps + i]
end
function Base.setindex!(x::TwoWaveguideTimestepView{T,false},input,i::Int) where T
    x.data[x.offset + (x.timeindex-1)*x.nsteps + i] = input
end
Base.size(x::TwoWaveguideTimestepView) = (x.nsteps,)
Base.eachindex(x::TwoWaveguideTimestepView) = 1:x.nsteps

"""
TwoWaveguideView{T} <: AbstractMatrix{T}

Structure for viewing state with one photon in Left and Right waveguide. 
"""
struct TwoWaveguideView{T,I} <: AbstractMatrix{T}
    data::AbstractArray{T}
    nsteps::Int
    offset::Int
end
TwoWaveguideView(data::AbstractArray{T},nsteps::Int,offset::Int,indexing) where {T} = TwoWaveguideView{T,indexing}(data,nsteps,offset)
Base.size(x::TwoWaveguideView) = (x.nsteps,x.nsteps)
function Base.getindex(x::TwoWaveguideView{T,true},i::Int,j::Int) where {T}
    x.data[x.offset + (i-1)*x.nsteps + j]
end
function Base.setindex!(x::TwoWaveguideView{T,true},input,i::Int,j::Int) where {T}
    x.data[x.offset + (i-1)*x.nsteps + j] = input
end
function Base.getindex(x::TwoWaveguideView{T,false},i::Int,j::Int) where {T}
    x.data[x.offset + (j-1)*x.nsteps + i]
end
function Base.setindex!(x::TwoWaveguideView{T,false},input,i::Int,j::Int) where {T}
    x.data[x.offset + (j-1)*x.nsteps + i] = input
end
