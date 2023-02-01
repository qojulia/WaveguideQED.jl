#Custom basis for Waveguide
mutable struct WaveguideBasis{P} <: QuantumOptics.Basis
    shape::Vector{Int}
    N::Int
    offset::Int
    nsteps::Int
    timeindex::Int
    dt::Float64
    function WaveguideBasis(N,times)
        dim = 0
        for i in 1:N
            dim = dim + length(times)^i - (length(times)^i-length(times))/2
        end
        new{N}([dim+1], dim, 0,length(times),1,times[2]-times[1])
    end
end

struct TwophotonView{T}
    data::T
    timeindex::Int
    nsteps::Int
end

function Base.getindex(x::TwophotonView,i::Int)
    if i<x.timeindex
        x.data[twophoton_index(i,x.nsteps,x.timeindex)]
    #XXX Double Check!
    elseif i == x.timeindex
        sqrt(2)*x.data[twophoton_index(x.timeindex,x.nsteps,x.timeindex)]
    else
        x.data[twophoton_index(x.timeindex,x.nsteps,i)]
    end
end


Base.eachindex(x::TwophotonView) = 1:x.nsteps

function Base.setindex!(x::TwophotonView,input,i::Int)
    if i<x.timeindex
        x.data[twophoton_index(i,x.nsteps,x.timeindex)] = input
    #XXX Double Check!
    elseif i == x.timeindex
        x.data[twophoton_index(x.timeindex,x.nsteps,x.timeindex)]  = sqrt(2)*input
    else
        x.data[twophoton_index(x.timeindex,x.nsteps,i)] = input
    end
end

function twophoton_index(j,nsteps,timeindex)
    1+nsteps+(j-1)*nsteps-((j-2)*(j-1))÷2+timeindex-j+1
end

#Function for creating "fockstate" or basically a very long vector for our basis.
# Namechange will happen in future, but for now kept for compatibility
function QuantumOpticsBase.:fockstate(::Type{T}, b::WaveguideBasis, n::Integer) where T
    @assert b.offset <= n <= b.N
    basisstate(T, b, n+1-b.offset)
end

function onephoton(b::WaveguideBasis,ξ)
    state = Ket(b)
    view = view_onephoton(state)
    view .= ξ
    normalize!(state)
    return state
end

function onephoton(b::WaveguideBasis,ξ::Function,times,args...)
    state = Ket(b)
    view = view_onephoton(state)
    view .= ξ.(times,args...)
    normalize!(state)
    return state
end

function twophoton(b::WaveguideBasis,ξ::Matrix)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = view(state.data,2+nsteps:Int(1+nsteps+nsteps*(nsteps+1)/2))
    for i in 1:nsteps
        for j in i:nsteps
            viewed_data[(i-1)*nsteps-((i-2)*(i-1))÷2+j-i+1] = ξ[i,j]
        end
    end
    normalize!(state)
    return state
end

function twophoton(b::WaveguideBasis,ξ::Function,times,args...)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = view(state.data,2+nsteps:Int(1+nsteps+nsteps*(nsteps+1)/2))
    for i in 1:nsteps
        for j in i:nsteps
            viewed_data[(i-1)*nsteps-Int((i-2)*(i-1)/2)+j-i+1] = ξ(times[i],times[j],args...)
        end
    end
    normalize!(state)
    return state
end

function get_waveguide_location(b::WaveguideBasis)
    return 1
end

function get_waveguide_location(b::CompositeBasis)
    return findall(x->typeof(x)<:WaveguideBasis,b.bases)
end

function view_waveguide(ψ::Ket)
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    view_waveguide(ψ,index)
end

function view_waveguide(ψ::Ket,index)
    view(reshape(ψ.data,Tuple(ψ.basis.shape)),index...)
end

function view_onephoton(ψ::Ket)
    viewed_data = view_waveguide(ψ::Ket)
    nsteps = get_nsteps(ψ.basis)
    return  view(viewed_data,2:nsteps+1)
end

function view_onephoton(ψ::Ket,index)
    viewed_data = view_waveguide(ψ::Ket,index)
    nsteps = get_nsteps(ψ.basis)
    return  view(viewed_data,2:nsteps+1)
end

function view_twophoton(ψ::Ket)
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    view_twophoton(ψ::Ket,index)
end



function view_twophoton(ψ::Ket,index)
    viewed_data = view_waveguide(ψ::Ket,index)
    nsteps = get_nsteps(ψ.basis)
    viewed_data = view(viewed_data,2+nsteps:Int(1+nsteps+nsteps*(nsteps+1)/2))
    output = zeros(ComplexF64,(nsteps,nsteps))
    for i in 1:nsteps
        for j in i:nsteps
            output[j,i] = viewed_data[(i-1)*nsteps-Int((i-2)*(i-1)/2)+j-i+1]
        end
    end
    output = output + output' - Diagonal(output)
end



function get_nsteps(basis::WaveguideBasis)
    basis.nsteps
end

function get_nsteps(basis::Basis)
    0
end

function get_nsteps(basis::CompositeBasis)
    for b in basis.bases
        if get_nsteps(b) != 0
            return get_nsteps(b)
        end
    end
end


function get_waveguidetimeindex(basis::WaveguideBasis)
    basis.timeindex
end

function get_waveguidetimeindex(basis::Basis)
    0
end

function get_waveguidetimeindex(basis::CompositeBasis)
    for b in basis.bases
        if get_waveguidetimeindex(b) != 0
            return get_waveguidetimeindex(b)
        end
    end
end