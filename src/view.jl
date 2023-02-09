"""
    TwophotonTimestepView{T}

Structure for viewing slice along same times in twophoton state.
"""
struct TwophotonTimestepView{T}
    data::T
    timeindex::Int
    nsteps::Int
end
function Base.getindex(x::TwophotonTimestepView,i::Int)
    if i<x.timeindex
        x.data[twophoton_index(i,x.nsteps,x.timeindex)]
    #XXX Double Check!
    elseif i == x.timeindex
        sqrt(2)*x.data[twophoton_index(x.timeindex,x.nsteps,x.timeindex)]
    else
        x.data[twophoton_index(x.timeindex,x.nsteps,i)]
    end
end
Base.eachindex(x::TwophotonTimestepView) = 1:x.nsteps
function Base.setindex!(x::TwophotonTimestepView,input,i::Int)
    if i<x.timeindex
        x.data[twophoton_index(i,x.nsteps,x.timeindex)] = input
    #XXX Double Check!
    elseif i == x.timeindex
        x.data[twophoton_index(x.timeindex,x.nsteps,x.timeindex)]  = sqrt(2)*input
    else
        x.data[twophoton_index(x.timeindex,x.nsteps,i)] = input
    end
end

"""
    TwophotonView{T} <: AbstractMatrix{T}

Structure for viewing twophoton state as symmetric matrix (only upper triangluar part is stored in memory). Returned from [`view_twophoton`](@ref). 
"""
struct TwophotonView{T} <: AbstractMatrix{T}
    data::AbstractArray{T}
    nsteps::Int
end
Base.size(x::TwophotonView) = (x.nsteps,x.nsteps)
function Base.getindex(x::TwophotonView,i::Int,j::Int)
    if i<j
        1/sqrt(2)*x.data[twophoton_index(i,x.nsteps,j)]
    #XXX Double Check!
    elseif i == j
        x.data[twophoton_index(i,x.nsteps,j)]
    else
        1/sqrt(2)*x.data[twophoton_index(j,x.nsteps,i)]
    end
end
function Base.setindex!(x::TwophotonView,input,i::Int,j::Int)
    if i<j
        x.data[twophoton_index(i,x.nsteps,j)] = sqrt(2)*input
    #XXX Double Check!
    elseif i == j
        x.data[twophoton_index(i,x.nsteps,j)] = input
    else
        println("i>j")
        x.data[twophoton_index(j,x.nsteps,i)] = sqrt(2)*input
    end
end
function twophoton_index(j,nsteps,timeindex)
    1+nsteps+(j-1)*nsteps-((j-2)*(j-1))÷2+timeindex-j+1
end