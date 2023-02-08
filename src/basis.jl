#Custom basis for Waveguide

"""
    WaveguideBasis(N, times)

Basis for time binned Waveguide where `N` is the number of photons in the waveguide.
Currently restricted to either 1 or 2. Times is timeinterval over which the photon state should be binned.


"""
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

"""
    onephoton(b::WaveguideBasis,ξ::Function,times,args...,norm=True)
    onephoton(b::WaveguideBasis,ξvec;norm=true)

Create a onephoton wavepacket of the form ``W^†(ξ) |0⟩ = \\int_{t_0}^{t_{end}} dt  ξ(t) w^†(t) |0⟩``. Here ``t_0=0`` and ``t_{end}`` is determined by the waveguide Basis BW.
ξ is a function that takes 

"""

function onephoton(b::WaveguideBasis,ξ::Function,times,args...; norm=true)
    state = Ket(b)
    view = view_onephoton(state)
    view .= ξ.(times,args...)
    if norm
        normalize!(state)
    end
    return state
end

function onephoton(b::WaveguideBasis,ξvec;norm=true)
    state = Ket(b)
    view = view_onephoton(state)
    view .= ξvec
    if norm
        normalize!(state)
    end
    return state
end

function twophoton(b::WaveguideBasis,ξ::Matrix;norm=true)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwophotonView(state.data,nsteps)
    for i in 1:nsteps
        for j in i:nsteps
            viewed_data[i,j] = ξ[i,j]
        end
    end
    if norm
        normalize!(state)
    end
    return state
end

function twophoton(b::WaveguideBasis,ξ::Function,times,args...;norm=true)
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwophotonView(state.data,nsteps)
    for i in 1:nsteps
        for j in i:nsteps
            viewed_data[i,j] = ξ(times[i],times[j],args...)
        end
    end
    if norm
        normalize!(state)
    end
    return state
end

function get_waveguide_location(b::WaveguideBasis)
    return 1
end

function get_waveguide_location(b::CompositeBasis)
    return findall(x->typeof(x)<:WaveguideBasis || typeof(x)<:OneTimeBasis,b.bases)[1]
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
    TwophotonView(viewed_data,nsteps)
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