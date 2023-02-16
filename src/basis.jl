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

Base.:(==)(b1::WaveguideBasis,b2::WaveguideBasis) = (b1.N==b2.N && b1.offset==b2.offset && b1.nsteps==b2.nsteps && b1.timeindex==b2.timeindex && b1.dt==b2.dt)

"""
    zerophoton(bw::WaveguideBasis)

Create a waveguide vacuum state |0⟩
"""
function zerophoton(b::WaveguideBasis)
    state = Ket(b)
    state.data[1] = 1
    return state
end


"""
    onephoton(b::WaveguideBasis,ξ::Function,times,args...,norm=True)
    onephoton(b::WaveguideBasis,ξvec;norm=true)

Create a onephoton wavepacket of the form ``W^†(ξ) |0⟩ = \\int_{t_0}^{t_{end}} dt  ξ(t) w^†(t) |0⟩``. Here ``t_0=0`` and ``t_{end}`` is determined by [`WaveguideBasis`](@ref).
ξ is a function evaluated as `ξ.(times,args...)`.
ξvec is a vector of length: `b.nsteps`.
If `norm==true` the state is normalized through `normalize!`.

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


"""
    twophoton(b::WaveguideBasis,ξ::Function,times,args...,norm=True)
    twophoton(b::WaveguideBasis,ξvec::Matrix;norm=true)

Create a twophoton wavepacket of the form ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  ξ(t,t') w^†(t)w^†(t') |0⟩``. Here ``t_0=0`` and ``t_{end}`` is determined by [`WaveguideBasis`](@ref).
ξ is a function evaluated as `ξ(t1,t2,args...)`. ξvec is a matrix of dimension: `(b.nsteps,b.nsteps)`, where `ξvec[i,j] = ξ(times[i],times[j])`, where times is defined in [`WaveguideBasis`](@ref).

"""
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


"""
    view_waveguide(ψ::ket)
    view_waveguide(ψ::ket,index)

View the Waveguide state given a state ψ containing a WaveguideBasis by returning `view(reshape(ψ.data,Tuple(ψ.basis.shape)),index...)`. If no index is provided the ground state is returned.
The index provided should be of the form `[:,i,j]` where `(:)` is at the location of the WaveguideBasis and i and j are indeces of other basises. See example: 

```
times=0:0.1:10
bw = WaveguideBasis(2,times)
bc1 = FockBasis(2)
bc2 = FockBasis(2)
ψ_waveguide = onephoton(bw,x->1)
ψ_total = ψ_waveguide ⊗ fockstate(bc1,1) ⊗ fockstate(bc2,1)
ψ_view = view_waveguide(ψ_total)
ψ_view_index = view_waveguide(ψ_total,[:,1,1])
ψ_view==ψ_view_index
```
"""
function view_waveguide(ψ::Ket)
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    view_waveguide(ψ,index)
end
function view_waveguide(ψ::Ket,index)
    view(reshape(ψ.data,Tuple(ψ.basis.shape)),index...)
end


"""
    view_onephoton(ψ::Ket)
    view_onephoton(ψ::Ket,index)

Return a view of the onephoton mode ``ξ(t)`` given an input state containing a onephoton waveguide state: ``\\int_{t_0}^{t_{end}} dt  ξ(t) w^†(t) |0⟩``
If no index is provided the ground state is returned. Index should follow same form outlined in [`view_waveguide`](@ref).

TO DO: PERHAPS CHANGE NAME?

"""
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


"""
    view_twophoton(ψ::Ket)
    view_twophoton(ψ::Ket,index)

Return a view of the twophoton mode ``ξ(t_1,t_2)`` given an input state containing a twophoton waveguide state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  ξ(t,t') w^†(t)w^†(t') |0⟩``
If no index is provided the ground state is returned. Index should follow same form outlined in [`view_waveguide`](@ref).

TO DO: PERHAPS CHANGE NAME?

"""
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


"""
    get_waveguide_location(basis::WaveguideBasis)
    get_waveguide_location(basis::CompositeBasis)

Return index of [`WaveguideBasis`](@ref) location in Hilbert space of basis b. `Btotal = BW ⊗ BC` where BW is a [`WaveguideBasis`](@ref) and `BC` some other basis then `get_waveguide_location(Btotal)` returns 1. 
While `Btotal = BC ⊗ BW` with `get_waveguide_location(Btotal)` returns 2.

"""
function get_waveguide_location(basis::WaveguideBasis)
    return 1
end
function get_waveguide_location(basis::CompositeBasis)
    return findall(x->typeof(x)<:WaveguideBasis,basis.bases)[1]
end


"""
    get_nsteps(basis::WaveguideBasis)
    get_nsteps(basis::Basis)
    get_nsteps(basis::CompositeBasis)

Return nsteps of [`WaveguideBasis`](@ref) given either a [`WaveguideBasis`](@ref) or a `CompositeBasis` containing a [`WaveguideBasis`](@ref)

"""
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


"""
    get_waveguidetimeindex(basis::WaveguideBasis)
    get_waveguidetimeindex(basis::Basis)
    get_waveguidetimeindex(basis::CompositeBasis)

Return timeindex of [`WaveguideBasis`](@ref) given either a [`WaveguideBasis`](@ref) or a `CompositeBasis` containing a [`WaveguideBasis`](@ref)

"""
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


"""
    get_waveguide_basis(basis::CompositeBasis)
    
Returns [`WaveguideBasis`](@ref) from `CompositeBasis.bases`

"""
function get_waveguide_basis(basis::CompositeBasis)
    basis.bases[findall(x->typeof(x)<:WaveguideBasis,basis.bases)]
end
