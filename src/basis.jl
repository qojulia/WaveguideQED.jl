"""
    WaveguideBasis(Np,Nw, times)
    aveguideBasis(Np, times)

Basis for time binned Waveguide where `Np` is the number of photons in the waveguide and `Nw` the number of waveguides (default is 1).
. Currently number of photons is restricted to either 1 or 2. Times is timeinterval over which the photon state should be binned.
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
            #Number of combinations of waveguides ([1,2],[1,3],[2,3] and so on)
            combinations = Nw*(Nw-1)/2
            dim = 1 + Nw*N + Nw*(N*(N+1))÷2 + N^2*combinations
        else
            error("Currently no more than two simultanoues photons are allowed")
        end
        new{Np,Nw}([dim+1], dim, 0,N)
    end
end
function WaveguideBasis(Np::Int,times)
    WaveguideBasis(Np,1,times)
end

#Type unions to dispatch on.
const SingleWaveguideBasis{Np} = Union{CompositeBasis{Vector{Int64}, T},WaveguideBasis{Np,1}} where {T<:Tuple{Vararg{Union{FockBasis,WaveguideBasis{Np,1}}}},Np}
const SingleWaveguideKet = Ket{T, Vector{ComplexF64}} where T<:SingleWaveguideBasis
const MultipleWaveguideBasis{Np,Nw} = Union{CompositeBasis{Vector{Int64},T},WaveguideBasis{Np,Nw},WaveguideBasis{Np,Nw}} where {T<:Tuple{Vararg{Union{FockBasis,WaveguideBasis{Np,Nw}}}},Np,Nw}
const MultipleWaveguideKet = Ket{T, Vector{ComplexF64}} where T<:MultipleWaveguideBasis



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


"""
    onephoton(b::WaveguideBasis{T,1},ξ::Function,times,args...,norm=True) where {T}
    onephoton(b::WaveguideBasis{T,1},ξvec;norm=true) where {T}
    onephoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Function,times,args...; norm=true) where {T,Nw}
    onephoton(b::WaveguideBasis{T,Nw},i,ξvec;norm=true) where {T,Nw}

Create a onephoton wavepacket in waveguide of the form ``W^\\dagger(\\xi) |0 \\rangle = \\int_{t_0}^{t_{end}} dt  \\xi(t) w_{\\mathrm{i}}^\\dagger(t) |\\emptyset \\rangle``.

# Arguments
- `b::WaveguideBasis{T,Nw}`: the basis of the waveguides, where T is the number of photons in the waveguides and Nw is the number of waveguides.
- `i::Int` (optional): the index of the waveguide where the photon is created. If not provided, i=1 is assumed.
- `ξ::Function` or `ξ::AbstractArray`: the wavefunction of the state. Can be a function with structure ξ.(times,args...) or an `AbstractArray` of length `b.nsteps`.
- `times`: A vector or range of times where the wavefunction is evaluated, used only if ξ is a function.
- `args...`: additional arguments to be passed to ξ if it is a function.
- `norm::Bool=true`: normalize the resulting wavepacket.


# Returns
- [`Ket(b)`]: a ket with the wavefunction defined above.

If `b` only contains one waveguide, the output wavefunction will contain on excitation in the Waveguide (i=1). If `b` contains multiple waveguides and only one index is given (i or j), then i=j is assumed.
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
function onephoton(b::WaveguideBasis{T,1},ξ::AbstractArray;norm=true) where {T}
    state = Ket(b)
    view = OnePhotonView(state)
    view .= ξ
    if norm
        normalize!(state)
    end
    return state
end
function onephoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Function,times,args...; norm=true) where {T,Nw}
    state = Ket(b)
    view = OnePhotonView(state,i)
    view .= ξ.(times,args...)
    if norm
        normalize!(state)
    end
    return state
end
function onephoton(b::WaveguideBasis{T,Nw},i::Int,ξ::AbstractArray;norm=true) where {T,Nw}
    state = Ket(b)
    view = OnePhotonView(state,i)
    view .= ξ
    if norm
        normalize!(state)
    end
    return state
end
onephoton(b::WaveguideBasis{T,Nw},ξ::Function,times,args...;norm=true) where {T,Nw} = onephoton(b,1,ξ,times,args...,norm=norm)
onephoton(b::WaveguideBasis{T,Nw},ξ::AbstractArray;norm=true) where {T,Nw} = onephoton(b,1,ξ,norm=norm)


"""
    twophoton(b::WaveguideBasis{T,1},ξ::Function,times,args...,norm=True) where {T}
    twophoton(b::WaveguideBasis{T,1},ξvec::Matrix;norm=true) where {T}
    twophoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Function,times,args...;norm=true) where {T,Nw}
    twophoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Matrix;norm=true) where {T,Nw}
    twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Function,times,args...;norm=true) where {T,Nw}
    twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Matrix;norm=true) where {T,Nw}
    twophoton(b::WaveguideBasis{T,Nw},WI::I,ξ::Matrix;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}}
    twophoton(b::WaveguideBasis{T,Nw},WI::I,ξ::Function,times,args...;norm=true)

Create a twophoton wavepacket of the form ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.

# Arguments
- `b::WaveguideBasis{T,Nw}`: the basis of the waveguides, where T is the number of photons in the waveguides and Nw is the number of waveguides.
- `i::Int` (optional): the index of the waveguide where the first photon is created. If not provided, i=1 is assumed.
- `j::Int` (optional): the index of the waveguide where the second photon is created. If not provided, j=i is assumed.
- `ξ::Function` or `ξ::Matrix`: the wavefunction of the state. Can be a function with structure ξ(times[l],times[m],args...) or a matrix of dimension `(b.nsteps, b.nsteps)`.
- `times`: A vector or range of times where the wavefunction is evaluated, used only if ξ is a function.
- `args...`: additional arguments to be passed to ξ if it is a function.
- `norm::Bool=true`: normalize the resulting wavepacket.


# Returns
- [`Ket(b)`]: a ket with the wavefunction defined above.

If `b` only contains one waveguide, the output wavefunction will contain two excitations in the same waveguide (i=j=1). If `b` contains multiple waveguides and only one index is given (i or j), then i=j is assumed, and two excitations in the same waveguide are returned.
`i` and `j` can also be given as a tuple or vector `WI= (i,j)` or `WI= [i,j]`
"""
function twophoton(b::WaveguideBasis{T,1},ξ::Function,times,args...;norm=true) where {T}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state.data,nsteps,nsteps+1)
    for l in 1:nsteps
        for m in l:nsteps
            viewed_data[l,m] = ξ(times[l],times[m],args...)
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
    for l in 1:nsteps
        for m in l:nsteps
            viewed_data[l,m] = ξ[l,m]
        end
    end
    if norm
        normalize!(state)
    end
    return state
end
function twophoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Function,times,args...;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,i)
    for l in 1:nsteps
        for m in l:nsteps
            viewed_data[l,m] = ξ(times[l],times[m],args...)
        end
    end
    if norm
        normalize!(state)
    end
    return state
end
function twophoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Matrix;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,i)
    for l in 1:nsteps
        for m in l:nsteps
            viewed_data[l,m] = ξ[l,m]
        end
    end
    if norm
        normalize!(state)
    end
    return state
end
function twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Function,times,args...;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,i,j)
    for l in 1:nsteps
        for m in l:nsteps
            viewed_data[l,m] = ξ(times[l],times[m],args...)
        end
    end
    if norm
        normalize!(state)
    end
    return state
end
function twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Matrix;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,i,j)
    for l in 1:nsteps
        for m in l:nsteps
            viewed_data[l,m] = ξ[l,m]
        end
    end
    if norm
        normalize!(state)
    end
    return state
end
twophoton(b::WaveguideBasis{T,Nw},idx::I,ξ::Matrix;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = twophoton(b,idx[1],idx[2],ξ,norm=norm)
twophoton(b::WaveguideBasis{T,Nw},idx::I,ξ::Function,times,args...;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = twophoton(b,idx[1],idx[2],ξ,times,args...,norm=norm)
twophoton(b::WaveguideBasis{T,Nw},ξ::Function,times,args...;norm=true) where {T,Nw} = twophoton(b,1,1,ξ,times,args...,norm=norm)
twophoton(b::WaveguideBasis{T,Nw},ξ::Matrix;norm=true) where {T,Nw} = twophoton(b,1,1,ξ,norm=norm)

"""
    view_waveguide(ψ::ket)
    view_waveguide(ψ::ket,index)

View the Waveguide state given a state ψ containing a WaveguideBasis by returning `view(reshape(ψ.data,Tuple(ψ.basis.shape)),index...)`. If no index is provided the ground state is returned.
The index provided should be of the form `[:,i,j]` where `(:)` is at the location of the WaveguideBasis and i and j are indeces of other basises. See example: 

```jldoctest
julia> using WaveguideQED; #hide
julia> using QuantumOptics; #hide
julia> times=0:0.1:10;
julia> bw = WaveguideBasis(2,times);
julia> bc1 = FockBasis(2);
julia> bc2 = FockBasis(2);
julia> ψ_waveguide = onephoton(bw,x->1);
julia> ψ_total = ψ_waveguide ⊗ fockstate(bc1,1) ⊗ fockstate(bc2,1);
julia> ψ_view = view_waveguide(ψ_total);
julia> ψ_view_index = view_waveguide(ψ_total,[:,1,1]);
julia> ψ_view==ψ_view_index
true
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

Return index of [`WaveguideBasis`](@ref) location in Hilbert space of basis b. `Btotal = waveguidebasis ⊗ otherbasis` where waveguidebasis is a [`WaveguideBasis`](@ref) and `otherbasis` some other basis then `get_waveguide_location(Btotal)` returns 1. 
While `Btotal = otherbasis ⊗ waveguidebasis` with `get_waveguide_location(Btotal)` returns 2.

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
    get_number_of_waveguides(basis::WaveguideBasis)
    get_number_of_waveguides(basis::Basis)
    get_number_of_waveguides(basis::CompositeBasis)

Return number of waveguides Nw of [`WaveguideBasis{Np,Nw}`](@ref) given either a [`WaveguideBasis{Np,Nw}`](@ref) or a `CompositeBasis` containing a [`WaveguideBasis{Np,Nw}`](@ref)
"""
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


