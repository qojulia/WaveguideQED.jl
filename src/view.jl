"""
    TwoPhotonTimestepView{T}

Structure for viewing slice along same times in twophoton state. Used in [`mul!`](@ref). 
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



"""
TwoWaveguideTimestepView{T} <:AbstractVector{T}

Structure for viewing slice along same times in twophoton states in two waveguides. Used in [`mul!`](@ref). 
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
