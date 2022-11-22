using Parameters
using NumericalIntegration
using DifferentialEquations

@with_kw mutable struct parameters
    γ::Float64 = 1.0
    δ::Float64 = 0.0
    σ::Float64 = 1.0
    t0::Float64 = 25.0
    N_cut::Int64 = 1.0
    times = 0.0:0.1:50.0
    x3::Float64 = 0.0
    δb::Float64 = 0.0
    Λn::Float64 = 0.0
end

function solve_onephotonwaveguide(param::parameters,ξvec)
    @unpack γ,δ,σ,t0,N_cut,times = param
    b_cw = cwbasis(times,1)
    a_dw = adw(b_cw)
    w_da = wda(b_cw)
    n = number(b_cw)
    #n = creation(b_cw)*annihilation(b_cw)
    dt = times[2]-times[1]
    #Create Hamiltonian for cavity waveguide system and onephoton conveyerbelt
    H = δ*n + (im*sqrt(γ/dt))*(a_dw-w_da)
    #Shape of input photon
    #Generate waveguide packet
    ξin = zero(b_cw)
    ξin.data[2] = OnePhoton(complex(0.0),sqrt(dt)*ξvec)
    #Check input normalized
    i1 = integrate(times,ξvec.*conj(ξvec))
    print("Input integral: $i1\n")
    
    cwsolver!(H,ξin)
    ξin.data[2].ξ01 = ξin.data[2].ξ01/sqrt(dt) 
    return ξin
end


function solve_twophotonwaveguide(param::parameters,ξvec)
    @unpack γ,δ,σ,t0,N_cut,times =param
    b_cw = cwbasis(times,2)
    dt = times[2] - times[1]
    #n = creation(b_cw)*annihilation(b_cw)
    n = number(b_cw)
    a_dw = adw(b_cw)
    w_da = wda(b_cw)
    #Create Hamiltonian for cavity waveguide system and onephoton conveyerbelt
    H = δ*n + (im*sqrt(γ/dt))*(a_dw-w_da)
    #print(H)
    #Shape of input photon
    #Generate waveguide packet
    ξin = zero(b_cw)
    ξin.data[3] = TwoPhoton(complex(0.0),zeros(ComplexF64,length(times)),ξvec)
    #Check input normalized
    
    i1 = zeros(ComplexF64,length(times))
    for j in 1:length(times)
        i1[j] = integrate(times,ξin.data[3].ξ02[j,:].*conj.(ξin.data[3].ξ02[j,:])/dt^2)
    end
    i2 = integrate(times,i1)
    print("Input integral: $i2\n") 

    cwsolver!(H,ξin)
    #Check output normalized
    for j in 1:length(times)
        i1[j] = integrate(times,ξin.data[3].ξ02[j,:].*conj.(ξin.data[3].ξ02[j,:])/dt^2)
    end
    i2 = integrate(times,i1)
    print("Output integral: $i2 \n")

    ξin.data[3].ξ02 = ξin.data[3].ξ02+ξin.data[3].ξ02'-Diagonal(ξin.data[3].ξ02)
    #ξin.data[3].ξ02 = ξin.data[3].ξ02-Diagonal(ξin.data[3].ξ02)
    return ξin
end

function solve_twophotonwaveguide_kerr(param::parameters,ξvec)
    @unpack γ,δ,σ,t0,N_cut,times,x3 = param
    b_cw = cwbasis(times,2)
    dt = times[2] - times[1]
    #n = creation(b_cw)*annihilation(b_cw)
    n = number(b_cw)
    a_dw = adw(b_cw)
    w_da = wda(b_cw)
    #Create Hamiltonian for cavity waveguide system and onephoton conveyerbelt
    H = δ*n + (im*sqrt(γ/dt))*(a_dw-w_da) + x3/4*(n*n-n)
    #print(H)
    #Shape of input photon
    #Generate waveguide packet
    ξin = zero(b_cw)
    ξin.data[3] = TwoPhoton(complex(0.0),zeros(ComplexF64,length(times)),ξvec)
    #Check input normalized
    
    i1 = zeros(ComplexF64,length(times))
    for j in 1:length(times)
        i1[j] = integrate(times,ξin.data[3].ξ02[j,:].*conj.(ξin.data[3].ξ02[j,:])/dt^2)
    end
    i2 = integrate(times,i1)
    print("Input integral: $i2\n") 

    cwsolver!(H,ξin)
    #Check output normalized
    for j in 1:length(times)
        i1[j] = integrate(times,ξin.data[3].ξ02[j,:].*conj.(ξin.data[3].ξ02[j,:])/dt^2)
    end
    i2 = integrate(times,i1)
    print("Output integral: $i2 \n")

    ξin.data[3].ξ02 = ξin.data[3].ξ02+ξin.data[3].ξ02'-Diagonal(ξin.data[3].ξ02)
    #ξin.data[3].ξ02 = ξin.data[3].ξ02-Diagonal(ξin.data[3].ξ02)
    return ξin
end

function solve_onephotonwaveguide_two_modes(param::parameters,ξvec)
    @unpack γ,δ,σ,t0,N_cut,times,δb,Λn = param
    b_cw = cwbasis(times,1)
    b_b = FockBasis(2)
    dt = times[2] - times[1]
    #n = creation(b_cw)*annihilation(b_cw)
    n = number(b_cw)
    a_dw = adw(b_cw)
    w_da = wda(b_cw)
    #Create Hamiltonian for cavity waveguide system and onephoton conveyerbelt
    H = δ*n + (im*sqrt(γ/dt))*(a_dw-w_da) + x3/4*(n*n-n)
    #print(H)
    #Shape of input photon
    #Generate waveguide packet
    ξin = zero(b_cw)
    ξin.data[3] = TwoPhoton(complex(0.0),zeros(ComplexF64,length(times)),ξvec)
    #Check input normalized
    
    i1 = zeros(ComplexF64,length(times))
    for j in 1:length(times)
        i1[j] = integrate(times,ξin.data[3].ξ02[j,:].*conj.(ξin.data[3].ξ02[j,:])/dt^2)
    end
    i2 = integrate(times,i1)
    print("Input integral: $i2\n") 

    cwsolver!(H,ξin)
    #Check output normalized
    for j in 1:length(times)
        i1[j] = integrate(times,ξin.data[3].ξ02[j,:].*conj.(ξin.data[3].ξ02[j,:])/dt^2)
    end
    i2 = integrate(times,i1)
    print("Output integral: $i2 \n")

    ξin.data[3].ξ02 = ξin.data[3].ξ02+ξin.data[3].ξ02'-Diagonal(ξin.data[3].ξ02)
    #ξin.data[3].ξ02 = ξin.data[3].ξ02-Diagonal(ξin.data[3].ξ02)
    return ξin
end



function solve_differentialeq(param::parameters,ξfun)
    @unpack γ,δ,σ,t0,N_cut,times = param
    p = [δ,γ]
    f(ψ,p,t) = complex(-(im*p[1]+p[2]/2)*ψ + sqrt(p[2])*ξfun(t,σ,t0))
    prob = ODEProblem(f,sqrt(p[2])*ξfun(times[1],σ,t0),(times[1],times[end]),p)
    solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
end
