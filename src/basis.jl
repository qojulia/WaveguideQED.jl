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
            dim = dim +length(times)^i
        end
        new{N}([dim+1], dim, 0,length(times),1,times[2]-times[1])
    end
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
    view .= sqrt(b.dt)*ξ
    return state
end

function twophoton(b::WaveguideBasis,ξ)
    state = Ket(b)
    view = view_twophoton(state)
    view .= sqrt(2)*b.dt*tril(ξ) + b.dt*Diagonal(ξ) - sqrt(2)*b.dt*Diagonal(ξ)
    return state
end

function onephoton(b::WaveguideBasis,ξ::Function,times,args...)
    state = Ket(b)
    view = view_onephoton(state)
    view .= sqrt(b.dt)*ξ.(times,args...)
    return state
end

function twophoton(b::WaveguideBasis,ξ::Function,times,args...)
    state = Ket(b)
    view = view_twophoton(state)
    mode = hcat([[ξ(t1,t2,args...) for t1 in times] for t2 in times]...)
    view  .= sqrt(2)*b.dt*tril(mode) + b.dt*Diagonal(mode) - sqrt(2)*b.dt*Diagonal(mode)
    return state
end

function get_waveguide_location(b::WaveguideBasis)
    return 1
end

function get_waveguide_location(b::CompositeBasis)
    return findall(x->typeof(x)<:WaveguideBasis,b.bases)
end

function view_waveguide(ψ::Ket)
    view = reshape(ψ.data,Tuple(ψ.basis.shape))
    loc = get_waveguide_location(ψ.basis)
    index = to_indices(view,Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape)))
    return view[index...]
end

function view_waveguide(ψ::Ket,index)
    view = reshape(ψ.data,Tuple(ψ.basis.shape))
    return view[index...]
end

function view_twophoton(ψ::Ket)
    viewed_data = view_waveguide(ψ::Ket)
    nsteps = get_nsteps(ψ.basis)
    return  reshape(view(viewed_data,2+nsteps:1+nsteps+nsteps^2),(nsteps,nsteps))
end

function view_twophoton(ψ::Ket{B,T}) where {B<:WaveguideBasis,T}
    nsteps = get_nsteps(ψ.basis)
    return  reshape(view(ψ.data,2+nsteps:1+nsteps+nsteps^2),(nsteps,nsteps))
end

function view_onephoton(ψ::Ket)
    viewed_data = view_waveguide(ψ::Ket)
    nsteps = get_nsteps(ψ.basis)
    return  view(viewed_data,2:nsteps+1)
end

function view_onephoton(ψ::Ket{B,T}) where {B<:WaveguideBasis,T}
    nsteps = get_nsteps(ψ.basis)
    return  view(ψ.data,2:nsteps+1)
end

function view_twophoton(ψ::Ket,index)
    viewed_data = view_waveguide(ψ::Ket,index)
    nsteps = get_nsteps(ψ.basis)
    return  reshape(view(viewed_data,2+nsteps:1+nsteps+nsteps^2),(nsteps,nsteps))
end

function view_singlephoton(ψ::Ket,index)
    viewed_data = view_waveguide(ψ::Ket,index)
    nsteps = get_nsteps(ψ.basis)
    return  view(viewed_data,2:nsteps+1)
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