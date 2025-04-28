using Test
using WaveguideQED
using LinearAlgebra
using QuantumOptics
include("helper_functions.jl")

@testset "Waveguide annihilation" begin
    times = 0:1:10
    dt = times[2] - times[1]

    bc = FockBasis(2)
    bw = WaveguideBasis(2,times)
    a = destroy(bc)
    ad = create(bc);
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w

    ξfun(t1,t2) = 1
    input = twophoton(bw,ξfun,norm=true)
    psi = fockstate(bc,0) ⊗ input
    psi_out = copy(psi)

    mul!(psi_out,adw,psi,1.0,0.0)
    ψ_single=OnePhotonView(psi_out,[2,:])
    testvec = ones(length(times)) .* 1/sqrt(get_nsteps(psi.basis)*(get_nsteps(psi.basis))÷2)
    @test isapprox(ψ_single,testvec,rtol=0.01)
end

@testset "Waveguide creation" begin
    times = 0:1:10
    dt = times[2] - times[1]

    bw = WaveguideBasis(2,times)
    bc = FockBasis(2)
    a = destroy(bc)
    ad = create(bc);
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w

    ξfun(t1) = 1
    input = onephoton(bw,ξfun,norm=true)
    psi = fockstate(bc,1) ⊗  input
    psi_out = copy(psi)
    mul!(psi_out,wda,psi,1,0.0)
    testvec = ones(length(times)) .* 1/sqrt(get_nsteps(psi.basis))
    testvec[get_waveguidetimeindex(wda)] *= 2
    ψ_double=TwoPhotonTimestepView(view_waveguide(psi_out),get_waveguidetimeindex(wda),get_nsteps(psi.basis),get_nsteps(psi.basis)+1)
    @test isapprox([ψ_double[i] for i in eachindex(ψ_double)],testvec)
end

@testset "1 photon 2 photon the same" begin
    times = 0:1:10
    dt = times[2] - times[1]

    bw = WaveguideBasis(2,times)
    bc = FockBasis(2)
    a = destroy(bc)
    ad = create(bc);
    n = ad*a ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w

    ξfun(t1) = 1
    input = onephoton(bw,ξfun,norm=true)
    psi = fockstate(bc,0) ⊗  input
    wda_out_2 = copy(psi)
    adw_out_2 = copy(psi)
    QuantumOpticsBase.mul!(wda_out_2,wda,psi,1,0.0)
    QuantumOpticsBase.mul!(adw_out_2,wda,psi,1,0.0)

    bw = WaveguideBasis(1,times)
    bc = FockBasis(2)
    a = destroy(bc)
    ad = create(bc);
    n = ad*a ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);
    wda = a ⊗ wd
    adw = ad ⊗ w

    input = onephoton(bw,ξfun,norm=true)
    psi = fockstate(bc,0) ⊗  input
    adw_out_1 = copy(psi)
    wda_out_1 = copy(psi)
    mul!(wda_out_1,wda,psi,1,0.0)
    mul!(adw_out_1,adw,psi,1,0.0)

    @test isapprox(OnePhotonView(wda_out_1),OnePhotonView(wda_out_2))
    @test isapprox(OnePhotonView(adw_out_1),OnePhotonView(adw_out_2))
    @test isapprox(wda_out_1.data[1],wda_out_2.data[1])
    @test isapprox(adw_out_1.data[1],adw_out_2.data[1])
end

@testset "Cavity Waveguide Operator" begin
    times = 0:1:10
    
    left_1photon = prep_order(1,1,1,1,[1,2,3])
    middle_1photon = prep_order(1,1,1,1,[2,3,1])
    right_1photon = prep_order(1,1,1,1,[3,1,2])
    
    @test unity_test(left_1photon)
    @test unity_test(middle_1photon)
    @test unity_test(right_1photon)

    left_2photon = prep_order(2,1,1,1,[1,2,3])
    middle_2photon = prep_order(2,1,1,1,[2,3,1])
    right_2photon = prep_order(2,1,1,1,[3,1,2])
    
    @test unity_test(left_2photon)
    @test unity_test(middle_2photon)
    @test unity_test(right_2photon)


    """
    dt = times[2] -times[1]
    bw = WaveguideBasis(2,times)
    bc = FockBasis(1)
    be = FockBasis(1)
    a = destroy(bc)
    ad = create(bc);
    n_l = identityoperator(be) ⊗ (ad*a) ⊗ identityoperator(bw)
    w = destroy(bw)
    wd = create(bw);

    wda_l = identityoperator(be) ⊗ emission(bc,bw)
    adw_l = identityoperator(be) ⊗ absorption(bc,bw)
    H_l = param.δ*n_l+im*sqrt(param.γ/dt)*(adw_l-wda_l)+param.x3/4*n_l*n_l-param.x3/4*n_l


    n_f = identityoperator(bw) ⊗ (ad*a) ⊗ identityoperator(be)
    wda_f = emission(bw,bc) ⊗ identityoperator(be)
    adw_f = absorption(bw,bc) ⊗ identityoperator(be)
    H_f = param.δ*n_f+im*sqrt(param.γ/dt)*(adw_f-wda_f)+param.x3/4*n_f*n_f-param.x3/4*n_f


    ξfun(t,σ,t0) = complex(1/(σ*sqrt(2*pi))*exp(-2*log(2)*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
    ξvec = ξfun.(times,1,1) * transpose(ξfun.(times,1,1))
    #Define initial state
    #ψ_cw = onephoton(bw,ξfun,times,param.σ,param.t0)
    ψ_cw = twophoton(bw,ξvec)

    psi_l = fockstate(be,0) ⊗ fockstate(bc,0) ⊗ ψ_cw
    out_l_emission = copy(psi_l)
    out_l_wda = copy(psi_l)
    out_l_absorption = copy(psi_l)
    out_l_adw = copy(psi_l)


    psi_f = ψ_cw ⊗ fockstate(bc,0) ⊗ fockstate(be,0)
    out_f_emission = copy(psi_f)
    out_f_wda = copy(psi_f)
    out_f_absorption = copy(psi_f)
    out_f_adw = copy(psi_f)


    set_waveguidetimeindex!([wda_l,wda_f,adw_l,adw_f],100)

    mul!(out_l_emission,wda_l,psi_l,1.0,0.0)
    mul!(out_f_emission,wda_f,psi_f,1.0,0.0)
    mul!(out_l_absorption,adw_l,psi_l,1.0,0.0)
    mul!(out_f_absorption,adw_f,psi_f,1.0,0.0)

    wda_l = identityoperator(be) ⊗ (a ⊗ wd)
    adw_l = identityoperator(be) ⊗ (ad ⊗ w)
    H_l = param.δ*n_l+im*sqrt(param.γ/dt)*(adw_l-wda_l)+param.x3/4*n_l*n_l-param.x3/4*n_l

    n_f = identityoperator(bw) ⊗ (ad*a) ⊗ identityoperator(be)
    wda_f = wd ⊗ a ⊗ identityoperator(be)
    adw_f = w ⊗ ad ⊗ identityoperator(be)
    H_f = param.δ*n_f+im*sqrt(param.γ/dt)*(adw_f-wda_f)+param.x3/4*n_f*n_f-param.x3/4*n_f

    set_waveguidetimeindex!([wda_l,wda_f,adw_l,adw_f],100)

    mul!(out_l_wda,wda_l,psi_l,1.0,0.0)
    mul!(out_f_wda,wda_f,psi_f,1.0,0.0)
    mul!(out_l_adw,adw_l,psi_l,1.0,0.0)
    mul!(out_f_adw,adw_f,psi_f,1.0,0.0)

    @test isapprox(out_l_wda,out_l_emission)
    @test isapprox(out_l_adw,out_l_absorption)
    @test isapprox(out_f_wda,out_f_emission)
    @test isapprox(out_f_adw,out_f_absorption)
    """
end


@testset "Multiple Waveguide Operators" begin
    for i in 1:3   
        left_1photon = prep_order(2,3,i,i,[1,2,3])
        middle_1photon = prep_order(2,3,i,i,[2,3,1])
        right_1photon = prep_order(2,3,i,i,[3,1,2])
        @test test_multiplewaveguides(left_1photon[1],3,i,[1,2,3],left_1photon[4][2],left_1photon[3][2])
        @test test_multiplewaveguides(middle_1photon[1],3,i,[2,3,1],middle_1photon[4][2],middle_1photon[3][2])
        @test test_multiplewaveguides(right_1photon[1],3,i,[3,1,2],right_1photon[4][2],right_1photon[3][2])
    end
    
    """
    times = 0:1:10
    dt = times[2] - times[1]
    nsteps = length(times)

    bw = WaveguideBasis(2,2,times)
    bc = FockBasis(2)
    wda_L = emission(bc,bw,1)
    wda_R = emission(bc,bw,2)
    adw_L = absorption(bc,bw,1)
    adw_R = absorption(bc,bw,2)

    tidx=5
    set_waveguidetimeindex!(wda_L,tidx)
    set_waveguidetimeindex!(wda_R,tidx)
    set_waveguidetimeindex!(adw_L,tidx)
    set_waveguidetimeindex!(adw_R,tidx)


    ξfun(t1) = 1
    psi = fockstate(bc,1) ⊗ onephoton(bw,1,ξfun,times)
    tmp = copy(psi)

    testvec = ones(nsteps) .* 1/sqrt(get_nsteps(psi.basis))
    testvec[get_waveguidetimeindex(wda_L)] *= 2
    mul!(tmp,wda_L,psi,1,0)
    psi_view = view_waveguide(tmp,[1,:])
    two = TwoPhotonTimestepView(psi_view,tidx,nsteps,1+2*nsteps);
    @test isapprox(two,testvec)
    mul!(tmp,adw_L,tmp,1,0)
    one = OnePhotonView(tmp,1,[2,:])
    @test isapprox(one,testvec)


    testvec = ones(nsteps) .* 1/sqrt(get_nsteps(psi.basis))
    mul!(tmp,wda_R,psi,1,0)
    psi_view = view_waveguide(tmp,[1,:])
    two = TwoWaveguideTimestepView(psi_view,tidx,nsteps,1+2*nsteps+(nsteps)*(nsteps+1),false);
    @test isapprox(two,testvec)
    mul!(tmp,adw_R,tmp,1,0)
    one = OnePhotonView(tmp,1,[2,:])
    @test isapprox(one,testvec)


    mul!(tmp,adw_L,psi,1,0)
    psi_view = view_waveguide(tmp,[3,:])
    @test isapprox(psi_view[1],1/sqrt(nsteps/2))
    mul!(tmp,wda_L,tmp,1,0)
    one = OnePhotonView(tmp,[2,:],1)
    @test isapprox(one[tidx],2/sqrt(nsteps))


    mul!(tmp,adw_R,psi,1,0)
    psi_view = view_waveguide(tmp,[3,:])
    @test isapprox(psi_view[1],0)
    mul!(tmp,wda_R,tmp,1,0)
    one = OnePhotonView(tmp,[2,:],2)
    @test isapprox(one[tidx],0)

    psi = fockstate(bc,1) ⊗ onephoton(bw,2,ξfun,times)

    testvec = ones(nsteps) .* 1/sqrt(get_nsteps(psi.basis))
    testvec[get_waveguidetimeindex(wda_L)] *= 2
    mul!(tmp,wda_R,psi,1,0)
    psi_view = view_waveguide(tmp,[1,:])
    two = TwoPhotonTimestepView(psi_view,tidx,nsteps,1+2*nsteps+(nsteps*(nsteps+1))÷2);
    @test isapprox(two,testvec)
    mul!(tmp,adw_R,tmp,1,0)
    one = OnePhotonView(tmp,[2,:],2)
    @test isapprox(one,testvec)

    testvec = ones(nsteps) .* 1/sqrt(get_nsteps(psi.basis))
    mul!(tmp,wda_L,psi,1,0)
    psi_view = view_waveguide(tmp,[1,:])
    two = TwoWaveguideTimestepView(psi_view,tidx,nsteps,1+2*nsteps+(nsteps)*(nsteps+1),true);
    @test isapprox(two,testvec)
    mul!(tmp,adw_L,tmp,1,0)
    one = OnePhotonView(tmp,[2,:],2)
    @test isapprox(one,testvec)


    mul!(tmp,adw_R,psi,1,0)
    psi_view = view_waveguide(tmp,[3,:])
    @test isapprox(psi_view[1],1/sqrt(nsteps/2))
    mul!(tmp,wda_R,tmp,1,0)
    one = OnePhotonView(tmp,[2,:],2)
    @test isapprox(one[tidx],2/sqrt(nsteps))


    mul!(tmp,adw_L,psi,1,0)
    psi_view = view_waveguide(tmp,[3,:])
    @test isapprox(psi_view[1],0)
    mul!(tmp,wda_L,tmp,1,0)
    one = OnePhotonView(tmp,[2,:],1)
    @test isapprox(one[tidx],0)
    """
end

@testset "Multiple Waveguide Interaction" begin
    combinations_list = [[1,2],[1,3],[2,3],[2,1],[3,2],[3,1]]
    for c in combinations_list  
        left_1photon = prep_order(2,3,c[1],c[2],[1,2,3])
        middle_1photon = prep_order(2,3,c[2],c[1],[2,3,1])
        right_1photon = prep_order(2,3,c[1],c[2],[3,1,2])
        @test test_interaction(left_1photon[1],3,c[1],c[2],[1,2,3],left_1photon[2][2])
        @test test_interaction(middle_1photon[1],3,c[2],c[1],[2,3,1],middle_1photon[2][2])
        @test test_interaction(right_1photon[1],3,c[1],c[2],[3,1,2],right_1photon[2][2])
        @test test_interaction(left_1photon[1],3,c[1],c[2],[1,2,3],left_1photon[2][1])
        @test test_interaction(middle_1photon[1],3,c[2],c[1],[2,3,1],middle_1photon[2][1])
        @test test_interaction(right_1photon[1],3,c[1],c[2],[3,1,2],right_1photon[2][1])
    end
end

@testset "Waveguide Delay Interaction" begin
    times = 0:1:10
    dt = times[2] - times[1]

    for i in 1:2
        bw = WaveguideBasis(i,2,times)
        
        tau = 5

        w1 = destroy(bw,1)
        wd1 = create(bw,1);
        w1_delay = destroy(bw,1;delay=tau/dt)
        wd1_delay = create(bw,1;delay=tau/dt);

        w2 = destroy(bw,2)
        wd2 = create(bw,2);
        w2_delay = destroy(bw,2;delay=tau/dt)
        wd2_delay = create(bw,2;delay=tau/dt);
        
        ξfun(t1) = 1
        input = Ket(bw)
        input.data .= 1
        normalize!(input)
        temp1 = copy(input)
        temp2 = copy(input)
        
        
        operators = [w1, wd1, w1_delay, wd1_delay, w2, wd2,w2_delay, wd2_delay]
        
        for o1 in operators
            for o2 in operators
                c1  = o1*o2
                c2 = LazyProduct([o1,o2])
                
                QuantumOptics.mul!(temp1,c1,input)
                QuantumOptics.mul!(temp2,c2,input)

                @test isapprox(temp1.data,temp2.data)
            end
        end
    end
    
end


@testset "NLevelWaveguideOperator with multiple bases" begin
    # Common parameters
    times = 0:0.1:10
    dt = times[2] - times[1]
    γL, γR = 1.0, 1.0
    
    @testset "TwoBases (NLevel + Waveguide)" begin
        # Set up bases
        bw = WaveguideBasis(2, 2, times)
        be = NLevelBasis(3)
        
        # Set up operators
        wh, wdh = destroy(bw, 1), create(bw, 1)
        wv, wdv = destroy(bw, 2), create(bw, 2)
        sh, sdh = transition(be, 1, 3), transition(be, 3, 1)
        sv, sdv = transition(be, 2, 3), transition(be, 3, 2)
        
        # Create identities for easier tensor products
        id_e = identityoperator(be)
        id_w = identityoperator(bw)
        
        # Test operators in different arrangements
        combined_basis = be ⊗ bw
        
        # Case 1: NLevel ⊗ Waveguide
        H1 = sdh ⊗ wv
        original_op1 = LazyTensor(combined_basis, combined_basis, (1, 2), (sdh, wv))
        
        # Create test states
        psi_in = Ket(be) ⊗ Ket(bw)
        
        # Initialize with random data for better testing
        psi_in.data .= rand(ComplexF64, length(psi_in.data))
        
        # Test Case 1
        psi_out1 = copy(psi_in)
        psi_ref1 = copy(psi_in)
        QuantumOptics.mul!(psi_out1, H1, psi_in)
        QuantumOptics.mul!(psi_ref1, original_op1, psi_in)
        @test isapprox(psi_out1.data, psi_ref1.data)



        # Test operators in different arrangements
        combined_basis = bw ⊗ be
        
        # Case 1: NLevel ⊗ Waveguide
        H1 =  wv ⊗ sdh
        original_op1 = LazyTensor(combined_basis, combined_basis, (1, 2), (wv, sdh))
        
        # Create test states
        psi_in = Ket(bw) ⊗ Ket(be)
        
        # Initialize with random data for better testing
        psi_in.data .= rand(ComplexF64, length(psi_in.data))
        
        # Test Case 1
        psi_out1 = copy(psi_in)
        psi_ref1 = copy(psi_in)
        QuantumOptics.mul!(psi_out1, H1, psi_in)
        QuantumOptics.mul!(psi_ref1, original_op1, psi_in)
        @test isapprox(psi_out1.data, psi_ref1.data)
    end

    # Test 1: Three bases - NLevel, Waveguide, and Fock basis as filler
    @testset "Three bases (NLevel + Waveguide + Fock)" begin
        # Set up bases
        bw = WaveguideBasis(2, 2, times)
        be = NLevelBasis(3)
        bf = FockBasis(4)  # Filler basis
        
        # Set up operators
        wh, wdh = destroy(bw, 1), create(bw, 1)
        wv, wdv = destroy(bw, 2), create(bw, 2)
        sh, sdh = transition(be, 1, 3), transition(be, 3, 1)
        sv, sdv = transition(be, 2, 3), transition(be, 3, 2)
        a, ad = destroy(bf), create(bf)
        
        # Create identities for easier tensor products
        id_e = identityoperator(be)
        id_w = identityoperator(bw)
        id_f = identityoperator(bf)
        
        # Test operators in different arrangements
        # Case 1: (NLevel ⊗ Fock) ⊗ Waveguide
        combined_basis = be ⊗ bf ⊗ bw
        H1 = (im*sqrt(γL)*sdh ⊗ id_f) ⊗ wv
        original_op1 = im*sqrt(γL)*LazyTensor(combined_basis, combined_basis, (1, 3), (sdh, wv))
        
        # Case 2: NLevel ⊗ (Fock ⊗ Waveguide)
        H2 = sdh ⊗ (id_f ⊗ wv)
        original_op2 = LazyTensor(combined_basis, combined_basis, (1, 3), (sdh, wv))
        
        # Case 3: Fock ⊗ (NLevel ⊗ Waveguide)
        combined_basis = bf  ⊗ be  ⊗ bw
        H3 = id_f ⊗ (sdh ⊗ wv)
        original_op3 = LazyTensor(combined_basis, combined_basis, (2, 3), (sdh, wv))
        
        # Create test states
        psi_in = Ket(be) ⊗ Ket(bf) ⊗ Ket(bw)
        psi_in.data .= 1.0
        
        # Test Case 1
        psi_out1 = copy(psi_in)
        psi_ref1 = copy(psi_in)
        QuantumOptics.mul!(psi_out1, H1, psi_in)
        QuantumOptics.mul!(psi_ref1, original_op1, psi_in)
        @test isapprox(psi_out1.data, psi_ref1.data)
        
        # Test Case 2
        psi_out2 = copy(psi_in)
        psi_ref2 = copy(psi_in)
        QuantumOptics.mul!(psi_out2, H2, psi_in)
        QuantumOptics.mul!(psi_ref2, original_op2, psi_in)
        @test isapprox(psi_out2.data, psi_ref2.data)
        
        # Test Case 3
        psi_in = Ket(bf) ⊗ Ket(be) ⊗ Ket(bw)
        psi_in.data .= 1.0
        psi_out3 = copy(psi_in)
        psi_ref3 = copy(psi_in)
        QuantumOptics.mul!(psi_out3, H3, psi_in)
        QuantumOptics.mul!(psi_ref3, original_op3, psi_in)
        @test isapprox(psi_out3.data, psi_ref3.data)
    end
    
    # Test 2: Four bases - NLevel, Waveguide, and two filler bases
    @testset "Four bases (NLevel + Waveguide + 2 fillers)" begin
        # Set up bases
        bw = WaveguideBasis(2, 2, times)
        be = NLevelBasis(3)
        bf1 = FockBasis(2)  # First filler basis
        bf2 = NLevelBasis(2)  # Second filler basis
        
        # Set up operators
        wh, wdh = destroy(bw, 1), create(bw, 1)
        wv, wdv = destroy(bw, 2), create(bw, 2)
        sh, sdh = transition(be, 1, 3), transition(be, 3, 1)
        sv, sdv = transition(be, 2, 3), transition(be, 3, 2)
        a, ad = destroy(bf1), create(bf1)
        s12, s21 = transition(bf2, 1, 2), transition(bf2, 2, 1)
        
        # Create identities for easier tensor products
        id_e = identityoperator(be)
        id_w = identityoperator(bw)
        id_f1 = identityoperator(bf1)
        id_f2 = identityoperator(bf2)
        
        # Test several complex arrangements
        combined_basis = be ⊗ bf1 ⊗ bf2 ⊗ bw
        
        # Case 1: (NLevel ⊗ Fock ⊗ NLevel) ⊗ Waveguide
        H1 = (sdh ⊗ id_f1 ⊗ id_f2) ⊗ wv
        original_op1 = LazyTensor(combined_basis, combined_basis, (1, 4), (sdh, wv))
        
        # Create test states
        psi_in = Ket(be) ⊗ Ket(bf1) ⊗ Ket(bf2) ⊗ Ket(bw)
        psi_in.data .= 1.0
        
        # Test Case 1
        psi_out1 = copy(psi_in)
        psi_ref1 = copy(psi_in)
        QuantumOptics.mul!(psi_out1, H1, psi_in)
        QuantumOptics.mul!(psi_ref1, original_op1, psi_in)
        @test isapprox(psi_out1.data, psi_ref1.data)
        


        # Test several complex arrangements
        combined_basis = bf1 ⊗ be ⊗ bf2 ⊗ bw
        
        # Case 1: (NLevel ⊗ Fock ⊗ NLevel) ⊗ Waveguide
        H1 = (id_f1 ⊗ sdh ⊗ id_f2) ⊗ wv
        original_op1 = LazyTensor(combined_basis, combined_basis, (2, 4), (sdh, wv))
        
        # Create test states
        psi_in = Ket(bf1) ⊗ Ket(be) ⊗ Ket(bf2) ⊗ Ket(bw)
        psi_in.data .= 1.0
        
        # Test Case 1
        psi_out1 = copy(psi_in)
        psi_ref1 = copy(psi_in)
        QuantumOptics.mul!(psi_out1, H1, psi_in)
        QuantumOptics.mul!(psi_ref1, original_op1, psi_in)
        @test isapprox(psi_out1.data, psi_ref1.data)

    end
    
end


