#Custom basis for Waveguide
mutable struct WaveguideBasis{P} <: QuantumOptics.Basis
    shape::Vector{Int}
    N::Int
    offset::Int
    nsteps::Int
    timeindex::Int
    function WaveguideBasis(N,times)
        dim = 0
        for i in 1:N
            dim = dim +length(times)^i
        end
        new{N}([dim+1], dim, 0,length(times),1)
    end
end

#Function for creating "fockstate" or basically a very long vector for our basis.
# Namechange will happen in future, but for now kept for compatibility
function QuantumOpticsBase.:fockstate(::Type{T}, b::WaveguideBasis, n::Integer) where T
    @assert b.offset <= n <= b.N
    basisstate(T, b, n+1-b.offset)
end


#TODO: Implement view which can show other than groundstate (|0,1_i,1_j>)
#TODO: Current implementation only works if WaveguideBasis is last basis.
#TODO: Use indexing used in mul! function defined in operators.jl
function view_twophoton(ψ::Ket)
    indeces = tr_indeces(ψ.basis) 
    viewed_data = view(ψ.data,indeces[1])
    for idx in indeces[2:end]
        viewed_data = view(viewed_data,idx)
    end
    nsteps = get_nsteps(ψ.basis)
    viewed_data = reshape(view(viewed_data,2+nsteps:1+nsteps+nsteps^2),(nsteps,nsteps))
end

function view_singlephoton(ψ::Ket)
    indeces = tr_indeces(ψ.basis) 
    viewed_data = view(ψ.data,indeces[1])
    for idx in indeces[2:end]
        viewed_data = view(viewed_data,idx)
    end
    nsteps = get_nsteps(ψ.basis)
    view(viewed_data,2:nsteps+1)
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

function tr_indeces(basis::WaveguideBasis)
    return [1:1:length(basis)]
end

function get_stride(basis::WaveguideBasis)
    return 1
end
function get_stride(basis::Basis)
    return basis.shape[1]
end

function tr_indeces(basis::CompositeBasis)
    indeces = Array{StepRange{Int64, Int64}}(undef,length(basis.bases))
    Hsize = Int(length(basis))
    for (i,b) in enumerate(basis.bases)
        stride = get_stride(b)
        indeces[i] = 1:Int(stride):Hsize
        Hsize = Int(Hsize/stride)
    end
    return indeces
end
