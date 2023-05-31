function _precompile_()
    order_lists = [[1,2,3],[2,1,3],[2,3,1]]
    for order in order_lists
        apply_all_op(prep_order(1,1,1,1,order))
        apply_all_op(prep_order(2,1,1,1,order))
        apply_all_op(prep_order(1,3,2,1,order))
        apply_all_op(prep_order(1,3,1,3,order))
        apply_all_op(prep_order(2,3,2,1,order))
        apply_all_op(prep_order(2,3,1,3,order))
    end

    times = 0:1:10
    dt = times[2] - times[1]
    bw = WaveguideBasis(1,1,times)
    bc = FockBasis(1)
    wda = emission(bc,bw)
    adw = absorption(bc,bw)
    H = im*sqrt(1/dt)*(adw-wda)

    #Define input onephoton state shape
    ξfun(t,σ,t0) = complex(sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2))
    ξvec = ξfun.(times,1,5)

    ψ_cw = onephoton(bw,ξvec)
    ψ_cw = onephoton(bw,ξfun,1,5)
    psi = fockstate(bc,0) ⊗  ψ_cw 

    #Run solvers.
    ψ = waveguide_evolution(times, psi, H)
    ψ = waveguide_montecarlo(times, psi, H,[destroy(bc) ⊗ identityoperator(bw)])
end

function prep_order(Np::Int,Nw::Int,wint_destroy::Int,wint_create::Int,order::Vector)
    times = 0:1:10
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
    return [b,nw_prod,nw,adw_prod,adw,wda_prod,wda]
end

function apply_all_op(oplist)
    psi = Ket(oplist[1])
    psi.data .=1
    tmp = copy(psi)
    for op in oplist[2:end]
        mul!(tmp,op,psi,1,0)
    end
end

using PrecompileTools

@setup_workload begin
    @compile_workload begin
        _precompile_()
    end
end