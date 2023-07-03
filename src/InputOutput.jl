"""
    WaveguideTransform{B1,B2,Np,idx} <: WaveguideOperator{B1,B2}

Operator structure for transforming an output state ``\\ket{\\psi}_{\\mathrm{out}} = \\mathbf{C} \\ket{\\tilde{\\psi}}_{\\mathrm{out}}``
Np is used to dispatch one or two photon routine and idx denotes the index of the waveguide the operator is acting on.
See also [`effective_hamiltonian`](@ref).

Examples:

```
times = 0:0.1:10
bw = WaveguideBasis(1,2,times)
C = 1/sqrt(2)*[1.0 -im;-im 1]
C_transform = WaveguideTransform(bw,C)
psi_tilde = ket(bw)
psi = C_transform*psi_tilde
```
"""
mutable struct WaveguideTransform{B1,B2,Np,Nw} <: WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    C::Matrix{ComplexF64}
    timeindex::Int
    function WaveguideTransform(bw::WaveguideBasis{Np,Nw},C::Matrix{ComplexF64}) where {Np,Nw}
        @assert Nw == size(C,1) == size(C,2)
        new{WaveguideBasis{Np,Nw},WaveguideBasis{Np,Nw},Np,Nw}(bw,bw,1,C,0)
    end
end

#Transform
function waveguide_mul!(result,a::WaveguideTransform{B,B,1,Nw},b,alpha,beta) where {B,Nw}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    fac = alpha*a.factor
    idx_end = (a.timeindex == 0) ? a.basis_l.nsteps : a.timeindex
    _onephotontransform!(result,a.C,b,fac,Nw,idx_end,a.basis_l.nsteps)
    return
end

function _onephotontransform!(result,C,b,fac,Nw,idx_end,nsteps)
    for i in 1:Nw
        @inbounds result[2+(i-1)*nsteps+idx_end:1+(i)*nsteps] .+= fac .* b[2+(i-1)*nsteps+idx_end:1+(i)*nsteps]
        for j in 1:Nw
            @inbounds result[2+(i-1)*nsteps:1+(i-1)*nsteps+idx_end] .+= C[i,j]*fac .* b[2+(j-1)*nsteps:1+(j-1)*nsteps+idx_end]
        end
    end
end

#Transform
function waveguide_mul!(result,a::WaveguideTransform{B,B,2,Nw},b,alpha,beta) where {B,Nw}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    fac = alpha*a.factor
    idx_end = a.timeindex == 0 ? a.basis_l.nsteps : a.timeindex
    _onephotontransform!(result,a.C,b,fac,Nw,idx_end,a.basis_l.nsteps)
    
    return
end


function io_relations_from_V(V,G)
    C = (I - im/2*V)*inv(I + im/2*V)
    correction = inv(I+im/2*V)
    G_corrected = correction*G
    Σ = 1/2*transpose(G)*correction*G
    C,V,correction,G_corrected,Σ
end
function io_relations_from_C(C,G) 
    V = -2*im*inv(I+C)*(I-C)
    Cout,V,correction,G_corrected,Σ = io_relations_from_V(V,G)
    @assert isapprox(Cout,C)
    C,V,correction,G_corrected,Σ
end

function Hinteraction_fromC(bw::WaveguideBasis{Np,Nw},bs,C,G) where {Np,Nw}
    C,V,correction,G_corrected,Σ = io_relations_from_C(C,G)
    n = identityoperator(bw)⊗number(bs)
    H = n*imag(Σ)*get_dt(bw)
    for i in 1:Nw
        H += (G_corrected[i])*absorption(bs,bw,i) + conj(G_corrected[i])*emission(bs,bw,i) 
    end
    H
end

"""
    effective_hamiltonian(bw::WaveguideBasis{Np,Nw},bs,C,G) where {Np,Nw}

Return the effective Hamiltonian for a waveguide system with basis `bw` coupled to a local emitter/cavity with the basis `bc`.
`C` and `G` here determine the input output relations of the coupling such that:
`` \\frac{d}{d t} a=-i\\left[a, H_{\\mathrm{c}}\\right]-\\Sigma a+\\mathbf{k}^T \\mathbf{W}_{\\mathrm{in}}``, and
`` \\mathbf{W}_{\\text {out }}(t)=\\mathbf{C} \\mathbf{W}_{\\mathrm{in}}(t)+a(t)\\mathbf{d}``. 
Here ``\\tilde{\\mathbf{d}} =  \\tilde{\\mathbf{k}} = -i \\left(\\left(\\mathbf{I}+\\frac{i}{2} \\mathbf{V}\\right)^{-1} \\right) \\mathbf{\\Gamma}``, where ``\\mathbf{\\Gamma}`` is determined by the vector `G` with a length equal to the number of waveguide in `bw`.
``\\mathbf{C}`` is here given by `C` which is of dimensions (nw,nw) with nw being the number of waveguides in `bw`.

The resulting state ``\\ket{\\tilde{\\psi}}_{\\mathrm{out}}`` from simulating with this Hamiltonian needs to be transformed using
[`WaveguideTransform`](@ref) such that ``\\ket{\\psi}_{\\mathrm{out}} = \\mathbf{C} \\ket{\\tilde{\\psi}}_{\\mathrm{out}}`` to get the correct input-output relations.

Examples:

```
times = 0:0.1:10
dt = times[2] - times[1]

bw = WaveguideBasis(1,2,times)
be = FockBasis(1)


C = [0 -1.0im; -1.0im 0]
γ=1
G = [sqrt(γ/dt),sqrt(γ/dt)]

H = effective_hamiltonian(bw,be,C,G)

ξfun(t1,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ^2)

ψ_in = onephoton(bw,1,ξfun,1,5) ⊗ fockstate(be,0)
ψ_out_tilde = waveguide_evolution(times,ψ_in,H)


C_transform = WaveguideTransform(bw,C) ⊗ identityoperator(be)
ψ_out = C_transform*ψ_out_tilde
```

"""
function effective_hamiltonian(bw::WaveguideBasis{Np,Nw},bs,C,G) where {Np,Nw}
    C,V,correction,G_corrected,Σ = io_relations_from_C(C,G)
    n = identityoperator(bw) ⊗ number(bs)
    H = n*imag(Σ)*get_dt(bw)
    for i in 1:Nw
        H += (G_corrected[i])*absorption(bw,bs,i) + conj(G_corrected[i])*emission(bw,bs,i) 
    end
    H
end
function effective_hamiltonian(bs,bw::WaveguideBasis{Np,Nw},C,G) where {Np,Nw}
    C,V,correction,G_corrected,Σ = io_relations_from_C(C,G)
    n = number(bs) ⊗ identityoperator(bw)
    H = n*imag(Σ)*get_dt(bw)
    for i in 1:Nw
        H += (G_corrected[i])*absorption(bs,bw,i) + conj(G_corrected[i])*emission(bs,bw,i) 
    end
    H
end

function  HinteractionV(bw::WaveguideBasis{Np,Nw},bs,V,G) where {Np,Nw}
    C,V,correction,G_corrected,Σ = io_relations_from_V(V,G)
    n = identityoperator(bw)⊗number(bs)
    H = n*imag(Σ)*get_dt(bw)
    for i in 1:Nw
        H += (G_corrected[i])*absorption(bw,bs,i) + conj(G_corrected[i])*emission(bw,bs,i) 
    end
    H
end


"""
    fftket(psi::Ket,idx)

Return the FFT ``\\xi(\\omega)`` of a onephoton wavefunction ``\\xi(t)``. Idx determines which waveguide.
"""
function fftket(psi::Ket,idx)
    xi = zeros(ComplexF64,length(OnePhotonView(psi,1)))
    xi .= OnePhotonView(psi,idx)
    spec = fftshift(fft(xi))
    freq =  fftshift(fftfreq(length(xi),1/get_dt(psi.basis))) .* (2*pi)
    return freq,spec
end


"""
    WaveguideTransformEvolution{B1,B2,Np,idx} <: WaveguideOperator{B1,B2}

Operator structure for transforming.
Np is used to dispatch one or two photon routine and idx denotes the index of the waveguide the operator is acting on. 
"""
mutable struct WaveguideTransformEvolution{B1,B2,Np,Nw} <: WaveguideOperator{B1,B2}
    basis_l::B1
    basis_r::B2
    factor::ComplexF64
    C::Matrix{ComplexF64}
    timeindex1::Int
    timeindex2::Int
    function WaveguideTransformEvolution(bw::WaveguideBasis{Np,Nw},C::Matrix{ComplexF64}) where {Np,Nw}
        @assert Nw == size(C,1) == size(C,2)
        new{WaveguideBasis{Np,Nw},WaveguideBasis{Np,Nw},Np,Nw}(bw,bw,1,C,0,0)
    end
end

#Transform
function waveguide_mul!(result,a::WaveguideTransformEvolution{B,B,1,Nw},b,alpha,beta) where {B,Nw}
    if iszero(beta)
        fill!(result, beta)
    elseif !isone(beta)
        rmul!(result, beta)
    end
    fac = alpha*a.factor
    idx1 = a.timeindex1
    idx2 = a.timeindex2
    if idx1>0 && idx2>0
        _onephotontransform_onetime!(result,a.C,b,fac,Nw,idx1,idx2,a.basis_l.nsteps)
    end
    return
end

function _onephotontransform_onetime!(result,C,b,fac,Nw,idx1,idx2,nsteps)
    for i in 1:Nw 
        for j in 1:Nw
            result[(i-1)*nsteps+idx1+1] += C[i,j]*fac*b[(j-1)*nsteps+idx2+1]
            if !iszero(C[i,j])
                result[(j-1)*nsteps+idx2+1] = 0
            end
        end
    end
end