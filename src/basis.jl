"""
    WaveguideBasis(N, times)

Basis for time binned Waveguide where `N` is the number of photons in the waveguide.
Currently restricted to either 1 or 2. Times is timeinterval over which the photon state should be binned.
"""
mutable struct WaveguideBasis{Np,Nw} <: QuantumOptics.Basis
    shape::Vector{Int}
    N::Int
    offset::Int
    nsteps::Int
    function WaveguideBasis(Np::Int,Nw::Int,times)
        dim = 0
        N = length(times)
        if Np == 1
            dim = 1 + N*Nw 
        elseif Np == 2
            combinations = Nw*(Nw-1)/2
            dim = 1 + Nw*N + Nw*(N*(N+1))÷2 + N^2*combinations
        else
            error("Currently no more than two simultanoues photons are allowed")
        end
        new{Np,Nw}([dim+1], dim, 0,N)
    end
end

const SingleWaveguideBasis{Np} = Union{CompositeBasis{Vector{Int64}, T},WaveguideBasis{Np,1}} where {T<:Tuple{Vararg{Union{FockBasis,WaveguideBasis{Np,1}}}},Np}
const SingleWaveguideKet = Ket{T, Vector{ComplexF64}} where T<:SingleWaveguideBasis
const MultipleWaveguideBasis{Np,Nw} = Union{CompositeBasis{Vector{Int64},T},WaveguideBasis{Np,Nw},WaveguideBasis{Np,Nw}} where {T<:Tuple{Vararg{Union{FockBasis,WaveguideBasis{Np,Nw}}}},Np,Nw}
const MultipleWaveguideKet = Ket{T, Vector{ComplexF64}} where T<:MultipleWaveguideBasis

function WaveguideBasis(Np::Int,times)
    WaveguideBasis(Np,1,times)
end

Base.:(==)(b1::WaveguideBasis,b2::WaveguideBasis) = (b1.N==b2.N && b1.offset==b2.offset && b1.nsteps==b2.nsteps)

"""
    zerophoton(bw::WaveguideBasis)

Create a waveguide vacuum state |0⟩
"""
function zerophoton(b::WaveguideBasis)
    state = Ket(b)
    state.data[1] = 1
    return state
end

function Noccupation(b::WaveguideBasis,i)
    state = Ket(b)
    state.data[i] = 1
    return state
end



"""
    onephoton(b::WaveguideBasis,ξ::Function,times,args...,norm=True)
    onephoton(b::WaveguideBasis,ξvec;norm=true)

Create a onephoton wavepacket of the form ``W^\\dagger(\\xi) |0⟩ = \\int_{t_0}^{t_{end}} dt  \\xi(t) w^\\dagger(t) |0⟩``. Here ``t_0=0`` and ``t_{end}`` is determined by [`WaveguideBasis`](@ref).
ξ is a function evaluated as `ξ.(times,args...)`.
ξvec is a vector of length: `b.nsteps`.
If `norm==true` the state is normalized through `normalize!`.

"""
function onephoton(b::WaveguideBasis{T,1},ξ::Function,times,args...; norm=true) where {T}
    state = Ket(b)
    view = OnePhotonView(state)
    view .= ξ.(times,args...)
    if norm
        normalize!(state)
    end
    return state
end
function onephoton(b::WaveguideBasis{T,1},ξvec;norm=true) where {T}
    state = Ket(b)
    view = OnePhotonView(state)
    view .= ξvec
    if norm
        normalize!(state)
    end
    return state
end
function onephoton(b::WaveguideBasis{T,Nw},idx::Int,ξ::Function,times,args...; norm=true) where {T,Nw}
    state = Ket(b)
    view = OnePhotonView(state,idx)
    view .= ξ.(times,args...)
    if norm
        normalize!(state)
    end
    return state
end
function onephoton(b::WaveguideBasis{T,Nw},idx,ξvec;norm=true) where {T,Nw}
    state = Ket(b)
    view = OnePhotonView(state,idx)
    view .= ξvec
    if norm
        normalize!(state)
    end
    return state
end


"""
    twophoton(b::WaveguideBasis,ξ::Function,times,args...,norm=True)
    twophoton(b::WaveguideBasis,ξvec::Matrix;norm=true)

Create a twophoton wavepacket of the form ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w^\\dagger(t)w^\\dagger(t') |0⟩``. Here ``t_0=0`` and ``t_{end}`` is determined by [`WaveguideBasis`](@ref).
ξ is a function evaluated as `ξ(t1,t2,args...)`. ξvec is a matrix of dimension: `(b.nsteps,b.nsteps)`, where `ξvec[i,j] = ξ(times[i],times[j])`, where times is defined in [`WaveguideBasis`](@ref).

"""
function twophoton(b::WaveguideBasis{T,1},ξ::Function,times,args...;norm=true) where {T}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state.data,nsteps,nsteps+1)
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
function twophoton(b::WaveguideBasis{T,1},ξ::Matrix;norm=true) where {T}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state.data,nsteps,nsteps+1)
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
function twophoton(b::WaveguideBasis{T,Nw},idx::Int,ξ::Function,times,args...;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,idx)
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
function twophoton(b::WaveguideBasis{T,Nw},idx::Int,ξ::Matrix;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,idx)
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
function twophoton(b::WaveguideBasis{T,Nw},WI1::Int,WI2::Int,ξ::Function,times,args...;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,WI1,WI2)
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
function twophoton(b::WaveguideBasis{T,Nw},WI1::Int,WI2::Int,ξ::Matrix;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,WI1,WI2)
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
twophoton(b::WaveguideBasis{T,Nw},WI::I,ξ::Matrix;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = twophoton(b,WI[1],WI[2],ξ,norm=norm)
twophoton(b::WaveguideBasis{T,Nw},WI::I,ξ::Function,times,args...;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = twophoton(b,WI[1],WI[2],ξ,times,args...,norm=norm)


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
    get_waveguide_location(basis::WaveguideBasis)
    get_waveguide_location(basis::CompositeBasis)

Return index of [`WaveguideBasis`](@ref) location in Hilbert space of basis b. `Btotal = BW ⊗ BC` where BW is a [`WaveguideBasis`](@ref) and `BC` some other basis then `get_waveguide_location(Btotal)` returns 1. 
While `Btotal = BC ⊗ BW` with `get_waveguide_location(Btotal)` returns 2.

"""
function get_waveguide_location(basis::WaveguideBasis)
    return 1
end
function get_waveguide_location(basis::CompositeBasis)
    return findall(x-> typeof(x)<:WaveguideBasis,basis.bases)[1]
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

function get_number_of_waveguides(x::WaveguideBasis{Np,Nw}) where {Np,Nw}
    Nw
end
function get_number_of_waveguides(basis::Basis)
    0
end
function get_number_of_waveguides(basis::CompositeBasis)
    for b in basis.bases
        if get_number_of_waveguides(b) != 0
            return get_number_of_waveguides(b)
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


