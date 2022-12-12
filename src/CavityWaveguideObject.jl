using QuantumOptics

export view_twophoton,get_cwbasis,get_woper,get_wdoper

#FOLLOWING IS COMMENTED SINCE NOT FINISHED

function view_twophoton(data::Vector{ComplexF64},nsteps::Int)
    reshape(view(data,2+nsteps:1+nsteps+nsteps^2),(nsteps,nsteps))
end

function view_twophoton(ψ::Ket,n::Int)
    ψ.basis.shape
end

"""
function view_twophoton(ψ::Ket,N::Int,nsteps::Int)
    
end
"""

function get_cwbasis(times,N)
    FockBasis(length(times)+(N-1)*length(times)^2)
end

function get_woper(b::Basis,N::Int,nsteps::Int,timeindex::Int)
    if N==1
        out = sparse(fockstate(b,0))
        start = sparse(dagger(fockstate(b,timeindex)))
        #+identityoperator(b) - tensor(out,dagger(out))- tensor(dagger(start),start)
        w = sparse(tensor(out,start))
    else
        out = sparse(fockstate(b,0))
        start = sparse(dagger(fockstate(b,timeindex)))
        w = sparse(tensor(out,start))
        
        for j in 1:timeindex-1
            out = sparse(fockstate(b,j))
            start = sparse(dagger(fockstate(b,nsteps+timeindex+(j-1)*nsteps)))
            w =w+sparse(tensor(out,start))
            #w = w + tensor(out,start)
        end
       
        for j in timeindex+1:nsteps
            out = sparse(fockstate(b,j))
            start = sparse(dagger(fockstate(b,nsteps+(timeindex-1)*nsteps+j)))
            w = w+sparse(tensor(out,start))
            end
        
        out = sparse(fockstate(b,timeindex))
        start = sparse(dagger(fockstate(b,nsteps+(timeindex-1)*nsteps+timeindex)))
        w = w + sqrt(2)*sparse(tensor(out,start))
    end
    return w
end

function get_wdoper(b::Basis,N::Int,nsteps::Int,timeindex::Int)
    if N==1
        out = sparse(fockstate(b,timeindex))
        start = sparse(dagger(fockstate(b,0)))
        w = sparse(tensor(out,start))
    else
        out = sparse(fockstate(b,timeindex))
        start = sparse(dagger(fockstate(b,0)))
        w = sparse(tensor(out,start))
        
        for j in 1:timeindex-1
            out = sparse(fockstate(b,nsteps+timeindex+(j-1)*nsteps))
            start = sparse(dagger(fockstate(b,j)))
            w =w+sparse(tensor(out,start))
            #w = w + tensor(out,start)
        end
        
        for j in timeindex+1:nsteps
            start = sparse(dagger(fockstate(b,j)))
            out = sparse(fockstate(b,nsteps+(timeindex-1)*nsteps+j))
            w = w+sparse(tensor(out,start))
            #w = w + tensor(out,start)
        end
        start = sparse(dagger(fockstate(b,timeindex)))
        out = sparse((fockstate(b,nsteps+(timeindex-1)*nsteps+timeindex)))
        w = w + sparse(tensor(out,start))
        #w = w + tensor(out,start)
        
    end
    return w
end


function solve(psi,f,times)
    dt = times[2] - times[1] 
    k1 = copy(psi)
    k2 = copy(psi)
    k3 = copy(psi)
    k4 = copy(psi)
    for i in 1:length(times)
        timestep!(psi,k1,k2,k3,k4,f,i,dt)
    end
    return psi
end

function timestep!(psi,k1,k2,k3,k4,f,t,dt)
    H = -im*f(t,psi)
    QuantumOptics.mul!(k1,H,psi,dt/2,0)
    #print(k1)

    k2.data .= k2.data+k1.data
    #print(k2)
    #print(k1)
    #QuantumOptics.mul!(k2,H,k2,dt/2,0)
    k2 = H*k2*dt/2

    k3.data .= k3.data+k2.data
    #QuantumOptics.mul!(k3,H,k3,dt,0)
    k3 = H*k3*dt

    k4.data .= k4.data+k3.data
    #QuantumOptics.mul!(k4,H,k4,dt/6,0)
    k4 = H*k4*dt/6

    k1 = 1/3*k1
    k2 = 2/3*k2
    k3 = 1/3*k3
    
    psi.data .= psi.data + k1.data+k2.data+k3.data+k4.data
    k1.data .= psi.data
    k2.data .= psi.data
    k3.data .= psi.data
    k4.data .= psi.data
end


"""
struct CWBasis <: Basis
    shape::Vector{Int}
    timeindex::Int
    NCutoff::Int
    nsteps::Int
end

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

