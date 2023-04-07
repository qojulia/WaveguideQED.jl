using Parameters
using DifferentialEquations

@with_kw mutable struct parameters
    γ::Float64 = 1.0
    δ::Float64 = 0.0
    σ::Float64 = 1.0
    t0::Float64 = 5.0
    times = 0.0:0.1:20.0
    x3::Float64 = 0.0
    δb::Float64 = 0.0
    Λn::Float64 = 0.0
end

function solve_differentialeq(param::parameters,ξfun)
    @unpack γ,δ,σ,t0,times = param
    p = [δ,γ]
    f(ψ,p,t) = complex(-(im*p[1]+p[2]/2)*ψ + sqrt(p[2])*ξfun(t,σ,t0))
    prob = ODEProblem(f,sqrt(p[2])*ξfun(times[1],σ,t0),(times[1],times[end]),p)
    solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,saveat=times)
end


@with_kw mutable struct BarretKokParameters
    γ::Float64 = 0.01
    g1::Float64 = 0.3
    g2::Float64 = 0.3
    κ1::Float64 = 1
    κ2::Float64 = 1
    δ::Float64 = 0.0
    td::Float64 = 5.0
    twait::Float64 = 30.0
    trelax::Float64 = 30.0
    times = 0.0:0.5:30.0
    montecarlo_repitions::Int = 10
end

function prep_fast(param)
    param.times = 0:0.5:param.twait
    dt = param.times[2] - param.times[1]
    param.td = param.twait

    #Create operators for two photons interacting with cavity
    bc = FockBasis(1)
    be = NLevelBasis(3)
    bw = WaveguideBasis(1,param.times)
    b_a = bw ⊗ bc ⊗ be
    a = embed(b_a,2,destroy(bc))
    ad = embed(b_a,2,create(bc))
    e = embed(b_a,3,transition(be,1,3))
    ed = embed(b_a,3,transition(be,3,1))
    wda = emission(bw,bc) ⊗ identityoperator(be)
    adw = absorption(bw,bc) ⊗ identityoperator(be)

    H_a = im*sqrt(param.κ1/2/dt)*(adw-wda) + param.g1/2*(ad*e + a*ed) 
    H_b = im*sqrt(param.κ2/2/dt)*(adw-wda)  + param.g2/2*(ad*e + a*ed)

    jump_a = sqrt(param.γ) *identityoperator(bw) ⊗ identityoperator(bc) ⊗ transition(be,1,3)
    jump_b = sqrt(param.γ) *identityoperator(bw) ⊗ identityoperator(bc) ⊗ transition(be,1,3) 

    wa = LazyTensor(b_a,1,destroy(bw))
    wb = LazyTensor(b_a,1,destroy(bw))

    n_w = LazyProduct(dagger(wa),wa)
    
    proj_up = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,2)
    proj_down = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,1)
    p1 = LazyTensorKet(proj_up,proj_down)
    p2 = LazyTensorKet(proj_down,proj_up)
    projector_plus = (p1+p2)/sqrt(2)
    projector_minus = (p1-p2)/sqrt(2)
  
    psi_a = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ (1/sqrt(2)*(nlevelstate(be,2)+nlevelstate(be,3)))
    psi_b = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ (1/sqrt(2)*(nlevelstate(be,2)+nlevelstate(be,3)))
    


    return (;(Base.@locals)...)
end


function prep_order(Np::Int,Nw::Int,wint_destroy::Int,wint_create::Int,order::Vector)
    times = 0:0.1:10
    bw = WaveguideBasis(Np,Nw,times)
    bc = FockBasis(2)
    be = FockBasis(1)
    basis_list = [bw,bc,be]
    b = tensor([basis_list[i] for i in order]...)

    bw_idx = findfirst(x->x==1,order)
    bc_idx = findfirst(x->x==2,order)

    w = embed(b,bw_idx,destroy(bw,wint_destroy))
    wd = embed(b,bw_idx,create(bw,wint_create))

    a = embed(b,bc_idx,destroy(bc))
    ad = embed(b,bc_idx,create(bc))

    nw_prod = LazyProduct(wd,w)
   
    adw_prod = LazyProduct(ad,w)
    wda_prod = LazyProduct(wd,a)

    nw = wd*w
    adw = tensor([i == bw_idx || i == bc_idx ? (i == bw_idx ? destroy(bw,wint_destroy) : create(bc)) : identityoperator(be) for i in 1:3]...)
    wda = tensor([i == bw_idx || i == bc_idx ? (i == bw_idx ? create(bw,wint_create) : destroy(bc)) : identityoperator(be) for i in 1:3]...)
    return [b,[nw_prod,nw],[adw_prod,adw],[wda_prod,wda]]
end
prep_order(Np,order::Vector) = prep_order(Np,1,1,order)

function unity_test(x)
    psi = Ket(x[1])
    psi.data .= 1
    psi1 = copy(psi)
    psi2 = copy(psi)
    for pairs in x[2:end]
        mul!(psi1,pairs[1],psi,1.123,0)
        mul!(psi2,pairs[2],psi,1.123,0)
        if !isapprox(psi1,psi2)
            return false
        end
    end
    return true
end

function test_multiplewaveguides(b,Nw,idx,order,wda1,adw1)
    bw_idx = findfirst(x->x==1,order)
    bc_idx = findfirst(x->x==2,order)
    ground_idx = [i == bw_idx ? (:) : 1 for i in 1:3]
    first_idx = [i == bw_idx ? (:) : i==bc_idx ? 2 : 1 for i in 1:3]
    second_idx = [i == bw_idx ? (:) : i==bc_idx ? 3 : 1 for i in 1:3]

    tidx = 5
    set_waveguidetimeindex!(adw1,tidx)
    set_waveguidetimeindex!(wda1,tidx)

    times = 0:0.1:10
    nsteps = length(times)
    dt = times[2] - times[1]
    ξfun(t1) = 1
    psi = tensor([i == bw_idx ? onephoton(b.bases[i],idx,ξfun,times) : i==bc_idx ? fockstate(b.bases[i],1) : fockstate(b.bases[i],0) for i in 1:3]...)
    tmp = copy(psi)
    
    testvec1 = ones(nsteps) .* 1/sqrt(get_nsteps(psi.basis))   

    
    testvec2 = ones(nsteps) .* 1/sqrt(get_nsteps(psi.basis))
    testvec2[get_waveguidetimeindex(wda1)] *= 2    

    mul!(tmp,wda1,psi,1,0)

    psi_view = view_waveguide(tmp,ground_idx)
    two = TwoPhotonTimestepView(psi_view,tidx,nsteps,1+Nw*nsteps+(idx-1)*(nsteps*(nsteps+1))÷2);
    if !isapprox(two,testvec2)
        return false
    end

    mul!(tmp,adw1,tmp,1,0)
    one = OnePhotonView(tmp,first_idx,idx)
    if !isapprox(one,testvec2)
        return false
    end

    mul!(tmp,adw1,psi,1,0)
    psi_view = view_waveguide(tmp,second_idx)
    if !isapprox(psi_view[1],1/sqrt(nsteps/2))
        return false
    end

    psi_view = view_waveguide(tmp,ground_idx)

    for k in filter(x -> x != idx, 1:Nw)
        psi = tensor([i == bw_idx ? onephoton(b.bases[i],k,ξfun,times) : i==bc_idx ? fockstate(b.bases[i],1) : fockstate(b.bases[i],0) for i in 1:3]...)
        mul!(tmp,wda1,psi,1,0)
        i,j = min(k,idx),max(k,idx)
        index = (i-1)*Nw + j - (i*(i+1))÷2
        twophotonview = TwoWaveguideTimestepView(psi_view,tidx,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx)
        if !isapprox(twophotonview,testvec1)
            return false
        end
    end
    return true
end


function test_interaction(b,Nw,idx1,idx2,order,wd2w1)
    bw_idx = findfirst(x->x==1,order)
    bc_idx = findfirst(x->x==2,order)
    ground_idx = [i == bw_idx ? (:) : 1 for i in 1:3]
    first_idx = [i == bw_idx ? (:) : i==bc_idx ? 2 : 1 for i in 1:3]
    second_idx = [i == bw_idx ? (:) : i==bc_idx ? 3 : 1 for i in 1:3]

    tidx = 5
    set_waveguidetimeindex!(wd2w1,tidx)

    times = 0:0.1:10
    nsteps = length(times)
    dt = times[2] - times[1]
    ξfun(t1) = 1
    psi = tensor([i == bw_idx ? onephoton(b.bases[i],idx1,ξfun,times) : fockstate(b.bases[i],0) for i in 1:3]...)
    tmp = copy(psi)
    
    testvec2 = ones(nsteps) .* sqrt(2)/get_nsteps(psi.basis)

    mul!(tmp,wd2w1,psi,1,0)

    psi_view = view_waveguide(tmp,ground_idx)
    if !isapprox(1/sqrt(get_nsteps(psi.basis)),psi_view[1+(idx2-1)*nsteps+tidx])
        return false
    end

    ξfun(t1,t2) = 1
    psi = tensor([i == bw_idx ? twophoton(b.bases[i],idx1,ξfun,times) : fockstate(b.bases[i],0) for i in 1:3]...)
    tmp = copy(psi)
    
    mul!(tmp,wd2w1,psi,1,0)
    psi_view = view_waveguide(tmp,ground_idx)
    
    i,j = min(idx1,idx2),max(idx1,idx2)
    index = (i-1)*Nw + j - (i*(i+1))÷2
    twophotonview = TwoWaveguideTimestepView(psi_view,tidx,nsteps,1+Nw*nsteps+Nw*(nsteps*(nsteps+1))÷2+(index-1)*nsteps^2,i==idx2)
    if !isapprox(twophotonview,testvec2)
        return false
    end
    return true
end