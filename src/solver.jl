"""
    waveguide_evolution(tspan, psi0, H; fout)

Integrate time-dependent Schroedinger equation to evolve states or compute propagators.
# Arguments
* `tspan`: Vector specifying the points of time for which output should be displayed.
* `psi0`: Initial state vector can only be a ket.
* `H`: Operator containing a [`WaveguideOperator`](@ref) either through a LazySum or LazyTensor.
* `fout=nothing`: If given, this function `fout(t, psi)` is called every time step. Example: `fout(t,psi) = expect(A,psi)` will return the epectation value of A at everytimestep.
   If `fout =1` the state psi is returned for all timesteps in a vector.
   ATTENTION: The state `psi` is neither normalized nor permanent! It is still in use by the ode solver and therefore must not be changed.


# Returns
* if `fout=nothing` the output of the solver will be the state `ψ` at the last timestep. 
* if `fout` is given a tuple with the state `ψ` at the last timestep and the output of `fout` is given. If `fout` returns a tuple the tuple will be flattened.
* if `fout = 1` `ψ` at all timesteps is returned.

# Examples 

* `fout(t,psi) = (expect(A,psi),expect(B,psi))` will result in  a tuple (ψ, ⟨A(t)⟩,⟨B(t)⟩), where `⟨A(t)⟩` is a vector with the expectation value of `A` as a function of time.

"""
function waveguide_evolution(times,psi,H;fout=nothing,kwargs...)
    ops = get_waveguide_operators(H)
    dt = get_dt(H.basis_l)
    nsteps = get_nsteps(H.basis_l)
    #Note that dt is added to the last input time to ensure that  
    tend = times[end]+dt
    times_sim = 0:dt:tend
    
    isapprox(norm(psi),1,rtol=10^(-6)) || @warn "Initial waveguidestate is not normalized. Consider passing norm=true to the state generation function."

    function get_hamiltonian(time,psi) 
        tidx = round(Int,time/dt,RoundDown) + 1
        set_waveguidetimeindex!(ops,tidx)
        return H
    end
    function eval_last_element(time,psi)
        if time >= tend
            return psi
        else
            return 0
        end
    end
    if fout === nothing
        tout, ψ = timeevolution.schroedinger_dynamic(times_sim, psi, get_hamiltonian,fout=eval_last_element;d_discontinuities=times_sim,kwargs...)
        return ψ[end]
    elseif fout == 1
        tout, ψ = timeevolution.schroedinger_dynamic(times_sim, psi, get_hamiltonian;d_discontinuities=times_sim,kwargs...)
        return ψ 
    else
        function feval(time,psi)
            if time == tend
                return (psi,fout(time,psi)...)
            else
                return (0,fout(time,psi)...)
            end
        end
        tout, ψ = timeevolution.schroedinger_dynamic(times_sim, psi, get_hamiltonian,fout=feval;d_discontinuities=times_sim,kwargs...)
        return (ψ[end][1], [[ψ[i][j] for i in 1:length(times_sim)-1] for j in 2:length(ψ[1])]...)
    end
end
    
"""
    waveguide_montecarlo(times,psi,H,J;fout=nothing)

See documentation for [`waveguide_evolution`](@ref) on how to define `fout`. J should be a list of collapse operators following documentation of [`timeevolution.mcwf_dynamic`](https://docs.qojulia.org/api/#QuantumOptics.timeevolution.mcwf_dynamic). 
"""
function waveguide_montecarlo(times,psi,H,J;fout=nothing,kwargs...)
    ops = get_waveguide_operators(H)
    dt = times[2] - times[1]
    tend = times[end]
    Jdagger = dagger.(J)
    
    isapprox(norm(psi),1,rtol=10^(-6)) || @warn "Initial waveguidestate is not normalized. Consider passing norm=true to the state generation function."

    function get_hamiltonian(time,psi)
        set_waveguidetimeindex!(ops,round(Int,time/dt,RoundUp)+1)
        return (H,J,Jdagger)
    end
    function eval_last_element(time,psi)
        if time == tend
            return psi
        else
            return 0
        end
    end
    if fout === nothing
        tout, ψ = timeevolution.mcwf_dynamic(times, psi, get_hamiltonian;fout=eval_last_element,kwargs...)
        return ψ[end]
    else
        function feval(time,psi)
            if time == tend
                return (psi,fout(time,psi)...)
            else
                return (0,fout(time,psi)...)
            end
        end
        tout, ψ = timeevolution.mcwf_dynamic(times, psi, get_hamiltonian;fout=feval,kwargs...)
        return (ψ[end][1], [[ψ[i][j] for i in 1:length(times)] for j in 2:length(ψ[1])]...)
    end
end

"""
    fast_unitary(times_eval,psi,H;order=2,fout=nothing)

See documentation for [`waveguide_evolution`](@ref) on how to define `fout`. J should be a list of collapse operators following documentation of [`timeevolution.mcwf_dynamic`](https://docs.qojulia.org/api/#QuantumOptics.timeevolution.mcwf_dynamic). 
"""
function fast_unitary(times_eval,psi,H;order=2,fout=nothing)
    if fout !== nothing
        fout_type = QuantumOptics.pure_inference(fout)
        output_container = Array{fout_type}(undef,length(times_eval))
    end
    isapprox(norm(psi),1,rtol=10^(-6)) || @warn "Initial waveguidestate is not normalized. Consider passing norm=true to the state generation function."
    dt = get_dt(H.basis_l)
    U = generate_unitary(H,dt,order)
    nsteps = min(get_nsteps(H.basis_l),round(Int,times_eval[end]/dt)+1)
    out = copy(psi)
    tmp = copy(psi)
    if (times_eval[2]-times_eval[1]) < dt
        @warn "Timestep of evaluation points is smaller than photon time binning dt. Dafaulting to sampling at every dt instead."
    end
    savefreq = round(Int,(times_eval[2]-times_eval[1])/dt)
    if fout === nothing
        for i in 1:nsteps
            set_waveguidetimeindex!(U,i)
            mul!(out,U,tmp,1,1)
            tmp.data .= out.data
        end
        return out
    else
        for i in 1:nsteps
            set_waveguidetimeindex!(U,i)
            mul!(out,U,tmp,1,1)
            tmp.data .= out.data
            if i%savefreq == 1
                output_container[i÷savefreq + 1] = fout(i*dt,out)
            end
        end
        return out,0:dt/savefreq:times_eval[end],output_container
    end
end

function generate_unitary(H,dt,order)
    U = (-im*dt)*H
    for i in 2:order
        U += (-im*dt)^i/factorial(i)*H^i
    end
    U
end