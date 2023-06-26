"""
    WaveguideTransform{B1,B2,Np,idx} <: WaveguideOperator{B1,B2}

Operator structure for transforming.
Np is used to dispatch one or two photon routine and idx denotes the index of the waveguide the operator is acting on. 
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

function  Hinteraction_fromC(bw::WaveguideBasis{Np,Nw},bs,C,G) where {Np,Nw}
    C,V,correction,G_corrected,Σ = io_relations_from_C(C,G)
    n = identityoperator(bw)⊗number(bs)
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
function  HinteractionV(bs,bw::WaveguideBasis{Np,Nw},V,G) where {Np,Nw}
    C,V,correction,G_corrected,Σ = io_relations_from_V(V,G)
    n = number(bs)⊗identityoperator(bw)
    H = n*imag(Σ)*get_dt(bw)
    for i in 1:Nw
        H += (G_corrected[i])*absorption(bs,bw,i) + conj(G_corrected[i])*emission(bs,bw,i) 
    end
    H
end

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
    #println(idx1)
    #println(idx1)
    #result .= b
    if idx1>0 && idx2>0
        _onephotontransform_onetime!(result,a.C,b,fac,Nw,idx1,idx2,a.basis_l.nsteps)
    end
    return
end

function _onephotontransform_onetime!(result,C,b,fac,Nw,idx1,idx2,nsteps)
    for i in 1:Nw
        #result[(i-1)*nsteps+idx1+1] = 0
        #result[2+(i-1)*nsteps+idx1:1+(i)*nsteps] .+= fac .* b[2+(i-1)*nsteps+idx1:1+(i)*nsteps]
        #result[2+(i-1)*nsteps:1+(i-1)*nsteps+idx1] .+= fac .* b[2+(i-1)*nsteps:1+(i-1)*nsteps+idx1]
        
        for j in 1:Nw
            result[(i-1)*nsteps+idx1+1] += C[i,j]*fac*b[(j-1)*nsteps+idx2+1]
            if !iszero(C[i,j])
                #println("zero")
                result[(j-1)*nsteps+idx2+1] = 0
            end
        end
    end
end