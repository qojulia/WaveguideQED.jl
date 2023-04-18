"""
    WaveguideBasis(Np,Nw, times)
    WaveguideBasis(Np, times)

Basis for time binned Waveguide where `Np` is the number of photons in the waveguide and `Nw` the number of waveguides (default is 1).
Currently number of photons is restricted to either 1 or 2. Times is timeinterval over which the photon state should be binned.
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

Base.:(==)(b1::WaveguideBasis,b2::WaveguideBasis) = (b1.N==b2.N && b1.offset==b2.offset && b1.nsteps==b2.nsteps)

#Type unions to dispatch on.
const SingleWaveguideBasis{Np} = Union{CompositeBasis{Vector{Int64}, T},WaveguideBasis{Np,1}} where {T<:Tuple{Vararg{Union{FockBasis,WaveguideBasis{Np,1}}}},Np}
const SingleWaveguideKet = Ket{T, Vector{ComplexF64}} where T<:SingleWaveguideBasis
const MultipleWaveguideBasis{Np,Nw} = Union{CompositeBasis{Vector{Int64},T},WaveguideBasis{Np,Nw},WaveguideBasis{Np,Nw}} where {T<:Tuple{Vararg{Union{FockBasis,WaveguideBasis{Np,Nw}}}},Np,Nw}
const MultipleWaveguideKet = Ket{T, Vector{ComplexF64}} where T<:MultipleWaveguideBasis


"""
    view_waveguide(ψ::ket)
    view_waveguide(ψ::ket,index)

View the Waveguide state given a state ψ containing a WaveguideBasis by returning `view(reshape(ψ.data,Tuple(ψ.basis.shape)),index...)`. If no index is provided the ground state is returned.
The index provided should be of the form `[:,i,j]` where `(:)` is at the location of the WaveguideBasis and i and j are indeces of other basises. See example: 

```@example vw
using QuantumOptics;
times=0:1:10;
bw = WaveguideBasis(1,times);
bc = FockBasis(2);
ψ_waveguide = Ket(bw,ones(length(times)+1));
ψ_total = ψ_waveguide ⊗ fockstate(bc,0) ⊗ fockstate(bc,0);
ψ_view = view_waveguide(ψ_total);
ψ_view_index = view_waveguide(ψ_total,[:,1,1]);
ψ_view==ψ_view_index
```

```@example vw
ψ_total = ψ_waveguide ⊗ fockstate(bc,2) ⊗ fockstate(bc,1);
view_waveguide(ψ_total,[:,3,2]) == ones(length(times))
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
    zerophoton(bw::WaveguideBasis)

Create a waveguide vacuum state ``| \\emptyset \rangle``
"""
function zerophoton(b::WaveguideBasis)
    state = Ket(b)
    state.data[1] = 1
    return state
end

"""
    OnePhotonView{T} <: AbstractVector{T}

Structure for viewing onephoton excitations in waveguides of the form ``W^\\dagger(\\xi) |0 \\rangle = \\int_{t_0}^{t_{end}} dt  \\xi(t) w_{\\mathrm{i}}^\\dagger(t) |\\emptyset \\rangle`` where ``i`` is the index of the waveguide.
See [`onephoton`](@ref) on how to create onephoton wavepackets and [`view_waveguide`](@ref) on how to index when there are multiple systems.

# Examples 
```@example onephotonview
using QuantumOptics;
times = 0:1:10;
bw = WaveguideBasis(1,times);
ψ1 = onephoton(bw,x->1,times,norm=false)
OnePhotonView(ψ1) == ones(length(times))
```

```@example onephotonview
bc = FockBasis(2);
ψ1Cavity = fockstate(bc,2,norm=false) ⊗ ψ1;
OnePhotonView(ψCavity,[3,:]) == ones(length(times))
```

```@example onephotonview
bw =  WaveguideBasis(1,3,times);
ψ2 = onephoton(bw,2,x->1,times,norm=false)
OnePhotonView(ψ,2) == ones(length(times))
```

```@example onephotonview
ψ2Cavity = fockstate(bc,2) ⊗ ψ2;
OnePhotonView(ψ2Cavity,2,[3,:]) == ones(length(times))
```
"""
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

ψ contains only a single waveguide and no waveguide index is needed.
No index of the system is provided and groundstate is assumed. 
Thus returns `OnePhotonView(ψ,index)` with `index = [1,:,1,...]` with `:` at the location of the waveguide and 1 in every other position. 
"""
function OnePhotonView(ψ::T) where {T<:SingleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    OnePhotonView(ψ,index)
end

"""
    OnePhotonView(ψ::T,index::I) where {T<:SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}

ψ contains only a single waveguide and no waveguide index is needed.
`index` should follow syntax outlined in [`view_waveguide`](@ref).
"""
function OnePhotonView(ψ::T,index::I) where {T<:SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    OnePhotonView(viewed_data,nsteps,1)
end

"""
    OnePhotonView(ψ::T,WI::Int)  where {T<:MultipleWaveguideKet}

ψ contains only a multiple waveguides and waveguide index `WI` is needed.
No index of the system is provided and groundstate is assumed (see [`view_waveguide`](@ref)). 
"""
function OnePhotonView(ψ::T,WI::Int)  where {T<:MultipleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    OnePhotonView(ψ,WI,index)
end

"""
    OnePhotonView(ψ::T,WI::Int,index::I) where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}

ψ contains only a multiple waveguides and waveguide index `WI` is needed.
`index` should follow syntax outlined in [`view_waveguide`](@ref).
"""
function OnePhotonView(ψ::T,WI::Int,index::I) where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    @assert WI <= get_number_of_waveguides(ψ.basis)
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    OnePhotonView(viewed_data,nsteps,1+nsteps*(WI-1))
end
OnePhotonView(ψ::T,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}} = throw(ArgumentError("When considering multiple waveguides please provide also index for the waveguide"))



"""
Create a onephoton wavepacket in waveguide of the form ``W^\\dagger(\\xi) |0 \\rangle = \\int_{t_0}^{t_{end}} dt  \\xi(t) w_{\\mathrm{i}}^\\dagger(t) |\\emptyset \\rangle`` where ``i`` is the index of the waveguide and return it as a `Ket`.
See also [`WaveguideBasis`](@ref) and [`OnePhotonView`](@ref) for how to view the state. 

# Examples
```@example onewaveguide
times = 1:1:10;
bw = WaveguideBasis(1,times);
ψ = onephoton(bw,x->1,times,norm=false);
OnePhotonView(ψ) == ones(length(times))
```
```@example onewaveguide
vec = collect(1:1:10);
ψ = onephoton(bw,vec);
OnePhotonView(ψ) == vec
```
```@example onewaveguide
bw = WaveguideBasis(1,3,times);
ψ = onephoton(bw,2,x->1,times);
OnePhotonView(ψ,2) == ones(length(times))
```
"""
function onephoton end

"""
    onephoton(b::WaveguideBasis{T,1},ξ::Function,times,args...,norm=True) where {T}
    
* Since `b` only contains a single waveguide, the index of the waveguide is not needed. 
* `ξ` should be broadcastable as ξ.(times,args...).
* `times`: A vector or range of times where the wavefunction is evaluated.
* `args...`: additional arguments to be passed to ξ if it is a function.
* `norm::Bool=true`: If true normalize the resulting wavepacket.
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

"""
    onephoton(b::WaveguideBasis{T,1},ξ::AbstractArray;norm=true) where {T}

* Since `b` only contains a single waveguide, the index of the waveguide is not needed. 
* `ξ` should be `AbstractArray` with length(ξ) == length(b.times)
* `norm::Bool=true`: If true normalize the resulting wavepacket.
"""
function onephoton(b::WaveguideBasis{T,1},ξ::AbstractArray;norm=true) where {T}
    state = Ket(b)
    view = OnePhotonView(state)
    view .= ξ
    if norm
        normalize!(state)
    end
    return state
end

"""
    onephoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Function,times,args...; norm=true) where {T,Nw}
    
* Since `b` contains `Nw` waveguides, the index of the waveguide needed.
* `i` is the index of the waveguide at which the onephoton wavepacket is created in
* `ξ` should be broadcastable as ξ.(times,args...).
* `times`: A vector or range of times where the wavefunction is evaluated.
* `args...`: additional arguments to be passed to ξ if it is a function.
* `norm::Bool=true`: If true normalize the resulting wavepacket.
"""
function onephoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Function,times,args...; norm=true) where {T,Nw}
    state = Ket(b)
    view = OnePhotonView(state,i)
    view .= ξ.(times,args...)
    if norm
        normalize!(state)
    end
    return state
end

"""
    onephoton(b::WaveguideBasis{T,Nw},i,ξvec;norm=true) where {T,Nw}

* Since `b` contains a `Nw` waveguides, the index of the waveguide needed.
* `i` is the index of the waveguide at which the onephoton wavepacket is created in
* `ξ` should be `AbstractArray` with length(ξ) == length(b.times)
* `norm::Bool=true`: If true normalize the resulting wavepacket.    
"""
function onephoton(b::WaveguideBasis{T,Nw},i::Int,ξ::AbstractArray;norm=true) where {T,Nw}
    state = Ket(b)
    view = OnePhotonView(state,i)
    view .= ξ
    if norm
        normalize!(state)
    end
    return state
end
onephoton(b::WaveguideBasis{T,Nw},ξ::Function,times,args...;norm=true) where {T,Nw} = throw(ArgumentError("WaveguideBasis contains multiple waveguides. Please provide the index of the waveguide in which the excitation should be created."))
onephoton(b::WaveguideBasis{T,Nw},ξ::AbstractArray;norm=true) where {T,Nw} = throw(ArgumentError("WaveguideBasis contains multiple waveguides. Please provide the index of the waveguide in which the excitation should be created."))




"""
    TwophotonView{T} <: AbstractMatrix{T}

Structure for viewing twophoton state of the form ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle`` where ``i`` and ``j`` are the indeces of the waveguide as matrix.
See also [`twophoton`](@ref) and [`view_waveguide`](@ref).

# Examples
Basic viewing:
```@example twophotview
using LinearAlgebra; #hide
times = 1:1:10;
bw = WaveguideBasis(2,times);
ψ = twophoton(bw,(t1,t2)->1,times,norm=false);
ψview = TwoPhotonView(ψ);
ψview == ones((length(times),length(times)))
```

Viewing state combined with other basis:
```@example twophotview
using QuantumOptics;
bc = FockBasis(2);
ψcombined = fockstate(bc,2) ⊗ ψ;
ψview = TwoPhotonView(ψcombined,[3,:]);
ψview == ones((length(times),length(times)))
```

Viewing twophoton state in waveguide 2 with multiple waveguides
```@example twophotview
bw = WaveguideBasis(2,3,times)
ψ = twophoton(bw,2,(t1,t2)->1,times,norm=false);
ψview = TwoPhotonView(ψ,2);
ψview == ones((length(times),length(times)))
```

Viewing twophoton state in waveguide 2 with multiple waveguides combined with other basis:
```@example twophotview
ψcombined = fockstate(bc,2) ⊗ ψ;
ψview = TwoPhotonView(ψcombined,2,[3,:]);
ψview == ones((length(times),length(times)))
```

Viewing twophotons across waveguide 1 and 2
```@example twophotview
bw = WaveguideBasis(2,3,times)
ψ = twophoton(bw,1,2,(t1,t2)->1,times,norm=false);
ψview = TwoPhotonView(ψ,1,2);
ψview == ones((length(times),length(times)))
```

Viewing twophotons across waveguide 1 and 2 combined with other basis:
```@example twophotview
ψcombined = fockstate(bc,2) ⊗ ψ;
ψview = TwoPhotonView(ψcombined,1,2,[3,:]);
ψview == ones((length(times),length(times)))
```
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

State `ψ` only contains one waveguide and no index provided so groundstate is assumed. See [`view_waveguide`](@ref).

"""
function TwoPhotonView(ψ::T) where {T <: SingleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,index)
end


"""
    TwoPhotonView(ψ::T,index::I) where {T <: SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}

State `ψ` only contains one waveguide. Index should follow syntax highlighted in [`view_waveguide`](@ref).
"""
function TwoPhotonView(ψ::T,index::I) where {T <: SingleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    TwoPhotonView(viewed_data,nsteps,nsteps+1)
end


"""
    TwoPhotonView(ψ::T,WI::Int) where {T <: MultipleWaveguideKet}

State `ψ` contains multiple waveguides and waveguide index `WI` required. 

Only one waveguide index `i=WI` means viewing the state ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle``.

No index provided so groundstate is assumed. See [`view_waveguide`](@ref).
"""
function TwoPhotonView(ψ::T,WI::Int) where {T <: MultipleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,WI,index)
end


"""
    TwoPhotonView(ψ::T,WI1::Int,WI2::Int) where {T <: MultipleWaveguideKet}

State `ψ` contains multiple waveguides.

Two waveguide indeces means viewing viewing the state ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle`` with `i = WI1` and `j = WI2`.

No index provided so groundstate is assumed. See [`view_waveguide`](@ref).
"""
function TwoPhotonView(ψ::T,WI1::Int,WI2::Int) where {T <: MultipleWaveguideKet}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,WI1,WI2,index)
end


"""
    TwoPhotonView(ψ::T,WI::I) where {T <: MultipleWaveguideKet,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}}

State `ψ` contains multiple waveguides.

Waveguide indeces provided as tuple or vector of length 2, means viewing viewing the state ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle`` with `i = WI[1]` and `j = WI[2]`.

No index provided so groundstate is assumed. See [`view_waveguide`](@ref).
"""
function TwoPhotonView(ψ::T,WI::I) where {T <: MultipleWaveguideKet,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}}
    loc = get_waveguide_location(ψ.basis)
    index = Tuple(i==loc[1] ? (:) : 1 for i in 1:length(ψ.basis.shape))
    TwoPhotonView(ψ,WI,index)
end


"""
    TwoPhotonView(ψ::T,WI::Int,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}

State `ψ` contains multiple waveguides.

Only one waveguide index `i=WI` means viewing the state ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle``.

Index should follow syntax highlighted in [`view_waveguide`](@ref).
"""
function TwoPhotonView(ψ::T,WI::Int,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    Nw = get_number_of_waveguides(ψ.basis)
    @assert WI<=Nw
    TwoPhotonView(viewed_data,nsteps,Nw*nsteps+1+(WI-1)*(nsteps*(nsteps+1)÷2))
end


"""
    TwoPhotonView(ψ::T,WI1::Int,WI2::Int,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}

State `ψ` contains multiple waveguides.

Two waveguide indeces means viewing viewing the state ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle`` with `i = WI1` and `j = WI2`.

Index should follow syntax highlighted in [`view_waveguide`](@ref).
"""
function TwoPhotonView(ψ::T,WI1::Int,WI2::Int,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}}}
    viewed_data = view_waveguide(ψ,index)
    nsteps = get_nsteps(ψ.basis)
    Nw = get_number_of_waveguides(ψ.basis)
    i,j = min(WI1,WI2),max(WI1,WI2)
    @assert i<=Nw && j<=Nw
    index = (i-1)*Nw + j - (i*(i+1))÷2
    TwoWaveguideView(viewed_data,nsteps,Nw*nsteps+1+(Nw)*(nsteps*(nsteps+1)÷2)+(index-1)*nsteps^2,i==WI1)
end


"""
    TwoPhotonView(ψ::T,WI::F,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}},F<:Union{Vector{Int64},Tuple{Vararg{Int64}}}}

State `ψ` contains multiple waveguides.

Waveguide indeces provided as tuple or vector of length 2, means viewing viewing the state ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle`` with `i = WI[1]` and `j = WI[2]`.

Index should follow syntax highlighted in [`view_waveguide`](@ref).

"""
TwoPhotonView(ψ::T,WI::F,index::I)  where {T<:MultipleWaveguideKet,I<:Union{Vector{Any},Vector{Int64},Tuple{Vararg{Union{Int64,Colon}}}},F<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = TwoPhotonView(ψ,WI[1],WI[2],index) 


"""
TwoWaveguideView{T} <: AbstractMatrix{T}

Structure for viewing state with one photon in waveguide i and j. Returned from [`TwoPhotonView`](@ref). See also [`TwoPhotonView`](@ref), [`twophoton`](@ref), and [`view_waveguide`](@ref)
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




"""
Create a twophoton wavepacket of the form ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle`` where ``i`` and ``j`` are the indeces of the waveguide and return it as a `Ket`.

See also [`WaveguideBasis`](@ref) and [`TwoPhotonView`](@ref) for how to view the state. 

#Examples 
Creating twophoton state with only one waveguide (using a function):
```@example twophotview
using LinearAlgebra; #hide
times = 1:1:10;
bw = WaveguideBasis(2,times);
ψ = twophoton(bw,(t1,t2)->1,times,norm=false);
ψview = TwoPhotonView(ψ);
ψview == ones((length(times),length(times)))

ψ = twophoton(bw,(t1,t2,arg1)->arg1,times,123,norm=false);
ψview = TwoPhotonView(ψ);
ψview == 123*ones((length(times),length(times)))

```

Creating twophoton state with only one waveguide (using a matrix):
```@example twophotview
ψ = twophoton(bw,ones((length(times),length(times))),norm=false);
ψview = TwoPhotonView(ψ);
ψview == ones((length(times),length(times)))
```

Creating twophoton state in waveguide 2 with multiple waveguides
```@example twophotview
bw = WaveguideBasis(2,3,times)
ψ = twophoton(bw,2,(t1,t2)->1,times,norm=false);
```

Creating twophoton state across waveguide 1 and 2
```@example twophotview
bw = WaveguideBasis(2,3,times)
ψ = twophoton(bw,1,2,(t1,t2)->1,times,norm=false);
```
"""
function twophoton end

"""
    twophoton(b::WaveguideBasis{T,1},ξ::Function,times,args...,norm=True) where {T}
 

# Arguments
- `b::WaveguideBasis{T,Nw}` contains only one waveguide and no index needed.
-  ξ given as a function should follow  `ξ(times[l],times[m],args...)`.
- `times`: A vector or range of length(times) = b.nsteps.
- `args...`: additional arguments to be passed to ξ if it is a function.
- `norm::Bool=true`: normalize the resulting wavepacket.
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

"""
    twophoton(b::WaveguideBasis{T,1},ξ::Matrix;norm=true) where {T}

# Arguments
- `b::WaveguideBasis{T,Nw}` contains only one waveguide and no index needed.
-  ξ given as a matrix of dimension `(b.nsteps, b.nsteps)`.
- `norm::Bool=true`: normalize the resulting wavepacket.
"""
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

"""
    twophoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Function,times,args...;norm=true) where {T,Nw}  

# Arguments
- `b::WaveguideBasis{T,Nw}` contains multiple waveguides and index i is need.
- `i` denotes the index of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle``.
-  ξ given as a function should follow  `ξ(times[l],times[m],args...)`.
- `times`: A vector or range of length(times) = b.nsteps.
- `args...`: additional arguments to be passed to ξ.
- `norm::Bool=true`: normalize the resulting wavepacket.
"""
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

"""
    twophoton(b::WaveguideBasis{T,Nw},i::Int,ξ::Matrix;norm=true) where {T,Nw}

# Arguments
- `b::WaveguideBasis{T,Nw}` contains multiple waveguides and index i is need.
- `i` denotes the index of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{i}}^\\dagger(t') |\\emptyset \\rangle``.
-  ξ given as a matrix of dimension `(b.nsteps, b.nsteps)`.
- `norm::Bool=true`: normalize the resulting wavepacket.
"""
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

"""
    twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Function,times,args...;norm=true) where {T,Nw}

# Arguments
- `b::WaveguideBasis{T,Nw}` contains multiple waveguides and index i is need.
- `i` denotes the index of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
- `j` denotes the index of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
-  ξ given as a function should follow  `ξ(times[l],times[m],args...)`.
- `times`: A vector or range of length(times) = b.nsteps.
- `args...`: additional arguments to be passed to ξ.
- `norm::Bool=true`: normalize the resulting wavepacket.
"""
function twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Function,times,args...;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,i,j)
    for l in 1:nsteps
        for m in 1:nsteps
            viewed_data[l,m] = ξ(times[l],times[m],args...)
        end
    end
    if norm
        normalize!(state)
    end
    return state
end

"""
    twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Matrix;norm=true) where {T,Nw}

# Arguments
- `b::WaveguideBasis{T,Nw}` contains multiple waveguides and index i is need.
- `i` denotes the index of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
- `j` denotes the index of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
-  ξ given as a matrix of dimension `(b.nsteps, b.nsteps)`.
- `norm::Bool=true`: normalize the resulting wavepacket.
"""
function twophoton(b::WaveguideBasis{T,Nw},i::Int,j::Int,ξ::Matrix;norm=true) where {T,Nw}
    state = Ket(b)
    nsteps = get_nsteps(b)
    viewed_data = TwoPhotonView(state,i,j)
    for l in 1:nsteps
        for m in 1:nsteps
            viewed_data[l,m] = ξ[l,m]
        end
    end
    if norm
        normalize!(state)
    end
    return state
end

"""
    twophoton(b::WaveguideBasis{T,Nw},idx::I,ξ::Matrix;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}}
# Arguments
- `b::WaveguideBasis{T,Nw}` contains multiple waveguides and index i is need.
- `I[1]` denotes the index ``i`` of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
- `I[2]` denotes the index ``j`` of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
-  ξ given as a matrix of dimension `(b.nsteps, b.nsteps)`.
- `norm::Bool=true`: normalize the resulting wavepacket.
"""
twophoton(b::WaveguideBasis{T,Nw},idx::I,ξ::Matrix;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = twophoton(b,idx[1],idx[2],ξ,norm=norm)

"""
    twophoton(b::WaveguideBasis{T,Nw},idx::I,ξ::Function,times,args...;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = twophoton(b,idx[1],idx[2],ξ,times,args...,norm=norm)

# Arguments
- `b::WaveguideBasis{T,Nw}` contains multiple waveguides and index i is need.
- `I[1]` denotes the index ``i`` of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
- `I[2]` denotes the index ``j`` of the waveguide in the twophoton state: ``\\int_{t_0}^{t_{end}} dt' \\int_{t_0}^{t_{end}} dt  \\xi(t,t') w_{\\mathrm{i}}^\\dagger(t)w_{\\mathrm{j}}^\\dagger(t') |\\emptyset \\rangle``.
-  ξ given as a function should follow  `ξ(times[l],times[m],args...)`.
- `times`: A vector or range of length(times) = b.nsteps.
- `args...`: additional arguments to be passed to ξ.
- `norm::Bool=true`: normalize the resulting wavepacket.
"""
twophoton(b::WaveguideBasis{T,Nw},idx::I,ξ::Function,times,args...;norm=true) where {T,Nw,I<:Union{Vector{Int64},Tuple{Vararg{Int64}}}} = twophoton(b,idx[1],idx[2],ξ,times,args...,norm=norm)

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


