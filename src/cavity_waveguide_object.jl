#################################################################################
#
# Cavity waveguide Module used for photon timebinning simulations
# Solver module for setting up parameters, and running the solver
#
# Matias Bundgaard-Nielsen
# DTU Fotonik 2022
#
#################################################################################
export generate_waveguide,waveguidestate,timestep!

import Base.+
import Base.-
import Base./
import Base.*

using QuantumOptics


abstract type waveguidestate end


#Structures containing the waveguidestate (1 or 2 photons)
mutable struct onephoton <:waveguidestate
    times
    N::Int
    ξ0::ComplexF64
    ξ1::Vector{ComplexF64}
    ψ::Ket
end

mutable struct twophoton <:waveguidestate
    times
    N::Int
    ξ0::ComplexF64
    ξ1::Vector{ComplexF64}
    ξ2::Matrix{ComplexF64}
    ψ::Ket
end

#Structures containing the differential of the waveguide state. 
#Requires less memory than a full waveguidestate because the derivative is only for one time
mutable struct onephoton_diff
    timeindex::Int
    ξ0::ComplexF64
    ξ1::ComplexF64    
end

mutable struct twophoton_diff
    timeindex::Int
    ξ0::ComplexF64
    ξ1::ComplexF64
    ξ2::Vector{ComplexF64}
end


#Function generatin waveguide
function generate_waveguide(
    times,
    N::Int,
    ψ0::Ket;
    ξ0::ComplexF64=complex(0.0),
    ξ1::Vector{ComplexF64}=zeros(ComplexF64,length(times)),
    ξ2::Matrix{ComplexF64}=zeros(ComplexF64,(length(times),length(times)))
    )

    if N==1
        onephoton(times,N,ξ0,ξ1,ψ0)
    elseif N==2
        twophoton(times,N,ξ0,ξ1,ξ2,ψ0)
    else
        print("Error: Currently no more than 2 photon states are supported")
    end
end

function generate_waveguide(state::waveguidestate)

    if state.N==1
        onephoton(state.times,state.N,complex(0.0),zeros(ComplexF64,length(state.times)),state.ψ)
    elseif state.N==2
        twophoton(state.times,state.N,complex(0.0),zeros(ComplexF64,length(state.times)),zeros(ComplexF64,(length(state.times),length(state.times))),state.ψ)
    else
        print("Error: Currently no more than 2 photon states are supported")
    end
end


#Operators that act on the cavity waveguidestate
function adw(state::onephoton,timeindex;g=1)
    return onephoton_diff(timeindex,complex(g*state.ξ1[timeindex]),complex(0.0))
end

function adw(state::twophoton,timeindex;g=1)
    diff_state = twophoton_diff(timeindex,complex(0.0),complex(0.0),zeros(ComplexF64,length(state.ξ2[timeindex,:])))
    diff_state.ξ0 =diff_state.ξ0 + sqrt(2)*complex(g*state.ξ1[timeindex])
    for j in 1:timeindex-1
        diff_state.ξ1 =diff_state.ξ1+complex(g*state.ξ2[timeindex,j])
    end
    for j in timeindex+1:length(state.ξ2[timeindex,:])
        diff_state.ξ1 =diff_state.ξ1+complex(g*state.ξ2[j,timeindex])
    end
    diff_state.ξ1 =diff_state.ξ1+complex(sqrt(2)*g*state.ξ2[timeindex,timeindex])
    return diff_state
end

function wda(state::onephoton,timeindex;g=1)
    return onephoton_diff(timeindex,complex(0.0),complex(g*state.ξ0))
end

function wda(state::twophoton,timeindex;g=1)
    diff_state = twophoton_diff(timeindex,complex(0.0),complex(0.0),zeros(ComplexF64,length(state.ξ2[timeindex,:])))
    diff_state.ξ1 = diff_state.ξ1 + sqrt(2)*complex(g*state.ξ0)
    for j in eachindex(state.ξ1)
        diff_state.ξ2[j] = diff_state.ξ2[j]+complex(g*state.ξ1[j])
    end
    diff_state.ξ2[timeindex] = diff_state.ξ2[timeindex]+(sqrt(2)-1)*complex(g*state.ξ1[timeindex])
    return diff_state
end

#Functions to get the current cavity state (old method). 
function get_state!(state::onephoton)
    state.ψ.data[2] = state.ξ0
    #When considering only a^dagger a, the occupation of the zero cavity state is unimportant
    c=1-state.ξ0*conj(state.ξ0)
    state.ψ.data[1] = c
    return c
end

function get_state!(state::twophoton)
    state.ψ.data[3] = state.ξ0
    if abs(state.ξ0)>1.0
    print("2-photon occ:  $(state.ξ0) \n")
    end
    c = sqrt(sum(state.ξ1.*conj(state.ξ1)))
    state.ψ.data[2] = c
    if abs(c)>1.0
    print("1-photon occ:  $(c) \n")
    end
    #When considering only a^dagger a, the occupation of the zero cavity state is unimportant
    state.ψ.data[1] = 1-c^2-state.ξ0^2
    return c
end

#Methods for adding waveguide differentials with waveguide differentials and states (used in timestep)
function +(state1::onephoton_diff,state2::onephoton_diff)
    return onephoton_diff(state1.timeindex,state1.ξ0+state2.ξ0,state1.ξ1+state2.ξ1)
end

function +(state1::twophoton_diff,state2::twophoton_diff)
    return twophoton_diff(state1.timeindex,state1.ξ0+state2.ξ0,state1.ξ1+state2.ξ1,state1.ξ2 .+ state2.ξ2 )
end


function *(state1::onephoton_diff,mul::ComplexF64)
    return onephoton_diff(state1.timeindex,state1.ξ0*mul,state1.ξ1*mul)
end
*(mul::ComplexF64,state1::onephoton_diff) = *(state1,mul) 
*(mul::T,state1::onephoton_diff) where T<:Number = *(state1,complex(float(mul))) 


function *(state1::twophoton_diff,mul::ComplexF64)
    return twophoton_diff(state1.timeindex,state1.ξ0*mul,state1.ξ1*mul,state1.ξ2 .* mul)
end
*(mul::ComplexF64,state1::twophoton_diff) = *(state1,mul) 
*(mul::T,state1::twophoton_diff) where T<:Number = *(state1,complex(float(mul))) 

/(state1::onephoton_diff,mul::ComplexF64,) = *(state1,1/mul) 
/(state1::twophoton_diff,mul::ComplexF64,) = *(state1,1/mul) 


function +(state1::onephoton,state_diff::onephoton_diff)
    state1.ξ0 = state1.ξ0 + state_diff.ξ0
    state1.ξ1[state_diff.timeindex] = state1.ξ1[state_diff.timeindex] + state_diff.ξ1
    return state1
end

function +(state1::twophoton,state_diff::twophoton_diff)
    state1.ξ0 = state1.ξ0 + state_diff.ξ0
    state1.ξ1[state_diff.timeindex] = state1.ξ1[state_diff.timeindex] + state_diff.ξ1 
    state1.ξ2[state_diff.timeindex,1:state_diff.timeindex] = state1.ξ2[state_diff.timeindex,1:state_diff.timeindex] + state_diff.ξ2[1:state_diff.timeindex] 
    state1.ξ2[state_diff.timeindex+1:end,state_diff.timeindex] = state1.ξ2[state_diff.timeindex+1:end,state_diff.timeindex]  + state_diff.ξ2[state_diff.timeindex+1:end] 
    return state1
end

function -(state1::onephoton,state_diff::onephoton_diff)
    state1.ξ0 = state1.ξ0 - state_diff.ξ0
    state1.ξ1[state_diff.timeindex] = state1.ξ1[state_diff.timeindex] - state_diff.ξ1 
    return state1
end

function -(state1::twophoton,state_diff::twophoton_diff)
    state1.ξ0 = state1.ξ0 - state_diff.ξ0
    state1.ξ1[state_diff.timeindex] = state1.ξ1[state_diff.timeindex] - state_diff.ξ1 
    state1.ξ2[state_diff.timeindex,1:state_diff.timeindex] = state1.ξ2[state_diff.timeindex,1:state_diff.timeindex] - state_diff.ξ2[1:state_diff.timeindex] 
    state1.ξ2[state_diff.timeindex+1:end,state_diff.timeindex] = state1.ξ2[state_diff.timeindex+1:end,state_diff.timeindex]  - state_diff.ξ2[state_diff.timeindex+1:end] 
    return state1
end


#Functions adding the cavity state back into the cavity waveguide state
function add_diff!(state::onephoton,c)
    state.ξ0 = state.ξ0 .+ state.ψ.data[2]
    state.ξ1 = state.ξ1 .+ state.ψ.data[1].*state.ξ1/c
end

function add_diff!(state::twophoton,c)
    state.ξ0 = state.ξ0 .+ state.ψ.data[3]
    if iszero(c)
        #print("c=0")
    else
        state.ξ1 = state.ξ1 .+ state.ψ.data[2]*state.ξ1/c
    end
    state.ξ2 = state.ξ2 .+ state.ψ.data[1]*state.ξ2
end


#The timestep algorithm using a RK4 method. Can be optimized by allocating objects for k1_w,k2_w,k3_w,k4_w
function timestep!(state,timeindex,γ,dt,H)
    # Cavity time evolution operator (first order)
    U = -im * H

    #Compute k1 derivatives
    c = get_state!(state)
    k1_w = adw(state,timeindex,g=sqrt(γ*dt)/6) + wda(state,timeindex,g=-sqrt(γ*dt)/6) 
    k1_c = U*state.ψ

    #Add k1 derivatives to cavity and waveguidestate to compute k2 state
    state.ψ = k1_c/2*dt
    add_diff!(state,c)  
    state = state + 3*k1_w
    
    c = get_state!(state)
    k2_w = adw(state,timeindex,g=sqrt(γ*dt)/3) + wda(state,timeindex,g=-sqrt(γ*dt)/3)
    k2_c = U*(state.ψ)
    
    #Substract k1 derivate and add k2
    state.ψ = k2_c/2*dt-k1_c/2*dt
    add_diff!(state,c)  
    state = state - 3*k1_w + 3/2*k2_w
    
    c = get_state!(state)
    k3_w = adw(state,timeindex,g=sqrt(γ*dt)/3) + wda(state,timeindex,g=-sqrt(γ*dt)/3)
    k3_c = U*(state.ψ)
    
    #Substract k2 derivate and add k3
    state.ψ=k3_c*dt-k2_c/2*dt
    add_diff!(state,c)  
    state = state - 3/2*k2_w + 3*k3_w
    
    c=get_state!(state)
    k4_w = adw(state,timeindex,g=sqrt(γ*dt)/6) + wda(state,timeindex,g=-sqrt(γ*dt)/6)
    k4_c = U*(state.ψ)

    #Substract k3 derivative
    state.ψ = -k3_c*dt
    add_diff!(state,c) 
    state = state -3*k3_w
    
    #Compute the total step using rk4
    state.ψ = 1/6*(k1_c+2*k2_c+2*k3_c+k4_c)*dt
    add_diff!(state,c)
    diff_state = k1_w+k2_w+k3_w+k4_w
    state+diff_state
end


export parameters,solve_differentialeq,solve_timebin,solve_onephotonwaveguide,solve_twophotonwaveguide

#Start of solver module
using QuantumOptics
using DifferentialEquations
using PyPlot
using NumericalIntegration
using LinearAlgebra

mutable struct parameters
    γ::Float64
    δ::Float64
    σ::Float64
    t0::Float64
    Np::Int64
    Nw::Int64
    N0::Int64
end



function timestep_solver!(H,ξin,γ)
    dt = ξin.times[2] - ξin.times[1] 
    ψ_out = Array{Ket}(undef,length(ξin.times))
    for i in 1:length(ξin.times)
        timestep!(ξin,i,γ,dt,H)
        ψ_out[i] = ξin.ψ
    end
    ξin.ξ1 = ξin.ξ1/sqrt(dt)
    i1 = integrate(ξin.times,ξin.ξ1.*conj.(ξin.ξ1))
    print("Integral new method: $i1 \n")
    return ψ_out
end



function solve_onephotonwaveguide(times,param::parameters) 
   
    dt =times[2] - times[1] 
    b_cavity = FockBasis(param.Np)
    a = destroy(b_cavity)
    ad = dagger(a)
    #Create Hamiltonian for cavity waveguide system and onephoton conveyerbelt
    H = param.δ*ad*a
    #Initial state of cavity:
    ψ0 = fockstate(b_cavity,param.N0)
    #Shape of input photon
    ξ(t::Float64,σ::Float64,t0::Float64) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
    #Generate waveguide packet
    ξin = generate_waveguide(times,param.Nw,ψ0,ξ1=ξ.(times,param.σ,param.t0)*sqrt(dt))
    
    #Check input normalized
    i1 = integrate(times,ξ.(times,param.σ,param.t0).^2)
    print("Input integral: $i1\n")

    #Call solver
    ψ_out = timestep_solver!(H,ξin,param.γ)
    return ξin,ψ_out
end

function solve_twophotonwaveguide(times,param::parameters) 
    dt =times[2] - times[1] 
    b_cavity = FockBasis(param.Np)
    a = destroy(b_cavity)
    ad = dagger(a)
    #Create Hamiltonian for cavity waveguide system and onephoton conveyerbelt
    H = param.δ*ad*a
    #Initial state of cavity:
    ψ0 = fockstate(b_cavity,param.N0)
    #Shape of input photon
    ξ(t::Float64,σ::Float64,t0::Float64) = complex(1/(σ*sqrt(2*pi))*exp(-1/2*(t-t0)^2/σ^2))/sqrt(0.2820947917738782)
    #Generate waveguide packet
    ξ2=tril(ξ.(times,param.σ,param.t0)*transpose(ξ.(times,param.σ,param.t0))*sqrt(2)*dt)
    #print(ξ2)
    ξin = generate_waveguide(times,param.Nw,ψ0,ξ2=ξ2)
    
    #i1 = integrate(times,ξ.(times,times,param.σ,param.t0).^2)
    #print("Input integral: $i1\n")

    ψ_out = timestep_solver!(H,ξin,param.γ)
    #ξin.ξ2 = ξin.ξ2[end:-1:1,:]
    #print(ξin.ξ2)
    ξin.ξ2 = ξin.ξ2+ξin.ξ2'-Diagonal(ξin.ξ2)
    #ξin.ξ2 = ξin.ξ2[end:-1:1,:]
    return ξin,ψ_out
end


function solve_timebin(ξ::Function,times,param::parameters) 
    dt = times[2] - times[1] 
    #Create FockBasis and operators for cavity containing N photons and waveguide with 1 photon
    b_cavity = FockBasis(param.Np)
    b_waveguide = FockBasis(1)
    a = destroy(b_cavity) ⊗ one(b_waveguide)
    ad = dagger(a)
    w = one(b_cavity) ⊗ destroy(b_waveguide)
    wd = dagger(w)
    b = b_cavity ⊗ b_waveguide
    #Create Hamiltonian for cavity waveguide system and onephoton conveyerbelt
    H = param.δ*ad*a + im*sqrt(param.γ/dt)*(ad*w-a*wd)
    #wg1 = generate_onephoton(ξ,times)
    U = - im * H*dt - 1/2*param.γ*dt*ad*a*w*wd
    #U = - im * H*dt - 1/2*H^2*dt^2
    #Initial state of cavity:
    ψ0 = 0*dm(fockstate(b_cavity,0)) ⊗ dm(fockstate(b_waveguide,1))
    #Output wavepacket
    output = Array{ComplexF64}(undef, length(times))
    #Loop over time and calculate derivative
    for i in 1:length(times)
        ψ0 = ψ0 + sqrt(dt)*ξ(times[i])*dm(fockstate(b_cavity,0)) ⊗ dm(fockstate(b_waveguide,1))
        diff = U*ψ0
        outcoupled_state = identityoperator(b_cavity) ⊗ dm(fockstate(b_waveguide,1))*(diff+ψ0)
        output[i] = tr(outcoupled_state)/sqrt(dt)
        ψ0 = ψ0 + diff - outcoupled_state
    end
    i1 = integrate(times,output .* conj.(output))
    print("Integral: $i1 \n")
    return output
end

function solve_differentialeq(ξ::Function, tspan, param::parameters)
    p = [param.δ,param.γ]
    f(ψ,p,t) = complex(-(im*p[1]+p[2]/2)*ψ + sqrt(p[2])*ξ(t,param.σ,param.t0))
    prob = ODEProblem(f,sqrt(p[2])*ξ(tspan[1],param.σ,param.t0),tspan,p)
    solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
end

