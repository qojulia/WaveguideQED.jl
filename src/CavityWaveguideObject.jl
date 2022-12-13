using QuantumOptics

export view_twophoton,get_cwbasis,get_woper,get_wdoper,WaveguideBasis,view_singlephoton

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


#ToDo implement view which can show other than groundstate (|0,1_i,1_j>)
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


function QuantumOpticsBase.:destroy(basis::WaveguideBasis{1})
    out = sparse(fockstate(basis,0))
    start = sparse(dagger(fockstate(basis,basis.timeindex)))
    sparse(tensor(out,start))
end

function QuantumOpticsBase.:destroy(basis::WaveguideBasis{2})
    out = sparse(fockstate(basis,0))
    start = sparse(dagger(fockstate(basis,basis.timeindex)))
    w = sparse(tensor(out,start))
    for j in 1:basis.timeindex-1
        out = sparse(fockstate(basis,j))
        start = sparse(dagger(fockstate(basis,basis.nsteps+basis.timeindex+(j-1)*basis.nsteps)))
        w =w+sparse(tensor(out,start))
    end
    for j in basis.timeindex+1:basis.nsteps
        out = sparse(fockstate(basis,j))
        start = sparse(dagger(fockstate(basis,basis.nsteps+(basis.timeindex-1)*basis.nsteps+j)))
        w = w+sparse(tensor(out,start))
    end
    out = sparse(fockstate(basis,basis.timeindex))
    start = sparse(dagger(fockstate(basis,basis.nsteps+(basis.timeindex-1)*basis.nsteps+basis.timeindex)))
    w = w + sqrt(2)*sparse(tensor(out,start))
end

function QuantumOpticsBase.:create(basis::WaveguideBasis)
    dagger(destroy(basis))
end

#FOLLOWING IS COMMENTED SINCE NOT FINISHED


"""


mutable struct CWState
    basis::CWBasis
    data::Vector{ComplexF64}
end

function copy(state::CWState)
    CWState(state.basis,copy(state.data))
end

struct CavWaveguideTensorState
    qo_b::Basis   
    cw_b::CWBasis
    data::Vector{ComplexF64}
end

struct CavWaveguideCompositeBasis <: Basis
    qo_b::Basis
    cw_b::CWBasis
end

tensor(qo_b::Basis,cw_b::CWBasis) = CavWaveguideCompositeBasis(qo_b,cw_b)

QuantumOptics.Base.:(⊗)(qo_b::Basis,cw_b::CWBasis) = tensor(qo_b,cw_b)

function view_twophoton(data::Vector{ComplexF64},b::CWBasis)
    return reshape(view(data,2+b.nsteps:b.nsteps^2),(b.nsteps,b.nsteps))
end

function w_zerophoton!(output::Vector{ComplexF64},input::Vector{ComplexF64},b::CWBasis;alpha=1,beta=1)
    "Nothing"
end

function w_onephoton!(output::Vector{ComplexF64},input::Vector{ComplexF64},b::CWBasis;alpha=1,beta=1)
    output[1] = beta*output[1] + complex(alpha*input[b.timeindex+1])
end

function w_twophoton!(output::Vector{ComplexF64},input::Vector{ComplexF64},b::CWBasis;alpha=1,beta=1)
    output[1] = beta*output[1] + complex(alpha*input[b.timeindex+1])
    two_photon_input = view_twophoton(input,b)
    for j in 1:b.timeindex-1
        output[j+1] = beta*output[j+1] + complex(alpha*two_photon_input[b.timeindex,j])
    end
    for j in b.timeindex+1:b.nsteps
        output[j+1] =beta*output[j+1]+complex(alpha*two_photon_input[j,b.timeindex])
    end
    output[b.timeindex+1] =beta*output[b.timeindex]+complex(sqrt(2)*alpha*two_photon_input[b.timeindex,b.timeindex])
end


function wd_zerophoton!(output::Vector{ComplexF64},input::Vector{ComplexF64},b::CWBasis;alpha=1,beta=1)
    "NOTHING"
end

function wd_onephoton!(output::Vector{ComplexF64},input::Vector{ComplexF64},b::CWBasis;alpha=1,beta=1)
    output[1+b.timeindex] = beta*output[1+b.timeindex] + complex(alpha*input[1])
end

function wd_twophoton!(output::Vector{ComplexF64},input::Vector{ComplexF64},b::CWBasis;alpha=1,beta=1)
    output[1+b.timeindex] = beta*output[1+b.timeindex] + complex(alpha*input[1])
    two_photon_output = view_twophoton(output,b)  
    for j in 1:b.timeindex-1
        two_photon_output[b.timeindex,j] = beta*two_photon_output[b.timeindex,j]+complex(alpha*input[j+1])
    end
    for j in b.timeindex+1:b.nsteps
        two_photon_output[j,b.timeindex] = beta*two_photon_output[j,b.timeindex] + complex(alpha*input[j+1])
    end
    two_photon_output[b.timeindex,b.timeindex] = beta*two_photon_output[b.timeindex,b.timeindex] + sqrt(2)*complex(alpha*input[b.timeindex+1])
end

"""

