    using WaveguideQED
    using QuantumOptics
    using PyPlot
    pygui(false)
    pygui(true)

    times = 0:0.1:10
    dt = times[2] - times[1]
    κ1 = 1
    κ2 = 1

    #Create operators for two    photons interacting with cavity
    be = FockBasis(1)
    bw = InputOutputWaveguideBasis(2,times)
    bw_single = InputOutputWaveguideBasis(1,times)
    b = bw ⊗ be
    b_single = bw_single ⊗ be
    a = embed(b,2,destroy(be))
    ad = embed(b,2,create(be))
    wdL = inputemission(bw,be)
    wL = inputabsorption(bw,be)
    wdR = outputemission(bw,be) 
    wR = outputabsorption(bw,be)

    a_single = embed(b_single,2,destroy(be))
    ad_single = embed(b_single,2,create(be))
    wdL_single = inputemission(bw_single,be)
    wL_single = inputabsorption(bw_single,be)
    wdR_single = outputemission(bw_single,be) 
    wR_single = outputabsorption(bw_single,be)

    H = im*sqrt(κ1/dt)*(wL-wdL) + im*sqrt(κ2/dt)*(wR-wdR)

    
    H_single = im*sqrt(κ1/dt)*(wL_single-wdL_single) + im*sqrt(κ2/dt)*(wR_single-wdR_single) 


    ξfun(t1,t2,σ1,σ2,t0) = sqrt(2/σ1)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ1^2)*sqrt(2/σ2)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t2-t0)^2/σ2^2)
    ξ_one_fun(t1,σ,t0) = sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t1-t0)^2/σ^2)

    psi_double_list = []
    psi_single_list = []
    widths = [0.8,1,1.5,2.9,8.3]/sqrt(4*log(2))
    t0 = 5

    for w in widths
        psi_in = twophoton(bw,:input,ξfun,times,w,w,t0) ⊗ fockstate(be,0)
        psi_in_single = onephoton(bw_single,:input,ξ_one_fun,times,w,t0) ⊗ fockstate(be,0)
        psi_out = waveguide_evolution(times,psi_in,H)
        psi_out_single = waveguide_evolution(times,psi_in_single,H_single)
        psi_R_scat = TwoPhotonView(psi_out,type=:input)
        push!(psi_double_list,psi_R_scat)

        psi_R_scat_single = zeros(ComplexF64,(length(times),length(times)))
        single_R = OnePhotonView(psi_out_single,type=:input)

        for i in eachindex(times)
            for j in eachindex(times)
                psi_R_scat_single[i,j] = single_R[i]*single_R[j]
            end
        end
        push!(psi_single_list,psi_R_scat_single)
    end

    fig,axs = subplots(2,5,figsize=(18,7))
    times = times .- 5
    xgrid = repeat(times',length(times),1)
    ygrid = repeat(times,1,length(times))
    for (i,ax) in enumerate(axs[2,:])
        cnt1= ax.contourf(xgrid,ygrid,psi_double_list[i].*conj(psi_double_list[i]),100,cmap="Blues")
        for c in cnt1.collections
            c.set_edgecolor("face")
        end
        ax.set_aspect("equal", "box")
        ax.set_xlabel(L"$t_1$")
        ax.set_ylabel(L"$t_2$")
        ax.set_title("")
    end

    for (i,ax) in enumerate(axs[1,:])
        cnt1= ax.contourf(xgrid,ygrid,psi_single_list[i].*conj(psi_single_list[i]),100,cmap="Reds")
        for c in cnt1.collections
            c.set_edgecolor("face")
        end
        ax.set_aspect("equal", "box")
        ax.set_yticks([])
        ax.set_xlabel(L"$t_1$")
    end
    plt.tight_layout()
    plt.savefig(pwd()*"/plots/lodahl_fig2.pdf")
