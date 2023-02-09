"""
    waveguide_evolution(tspan, psi0, H; fout)

Integrate time-dependent Schroedinger equation to evolve states or compute propagators.
# Arguments
* `tspan`: Vector specifying the points of time for which output should be displayed.
* `psi0`: Initial state vector can only be a ket.
* `H`: Operator containing a [`WaveguideOperator`](@ref) either through a LazySum or LazyTensor.
* `fout=nothing`: If given, this function `fout(t, psi)` is called every time step. Example: `fout(t,psi) = expect(A,psi)` will return the epectation value of A at everytimestep. 
   ATTENTION: The state `psi` is neither normalized nor permanent! It is still in use by the ode solver and therefore must not be changed.

# Output
* if `fout=nothing` the output of the solver will be the state `ψ` at the last timestep. 
* if `fout` is given a tuple with the state `ψ` at the last timestep and the output of `fout` is given. If `fout` returns a tuple the tuple will be flattened.
Example `fout(t,psi) = (expect(A,psi),expect(B,psi))` will result in  a tuple (ψ, ⟨A(t)⟩,⟨B(t)⟩), where `⟨A(t)⟩` is a vector with the expectation value of `A` as a function of time.
"""
function waveguide_evolution(times,psi,H;fout=nothing)
    basis = get_waveguide_basis(psi.basis)
    dt = times[2] - times[1]
    tend = times[end]
    function get_hamiltonian(time,psi)
        #index = findlast(times .<= time)
        #basis.timeindex = index
        basis.timeindex = round(Int,time/dt,RoundUp)+1
        return H
    end
    function eval_last_element(time,psi)
        if time == tend
            return psi
        else
            return 0
        end
    end
    if fout === nothing
        tout, ψ = timeevolution.schroedinger_dynamic(times, psi, get_hamiltonian,fout=eval_last_element)
        return ψ[end]
    else
        function feval(time,psi)
            if time == tend
                return (psi,fout(time,psi)...)
            else
                return (0,fout(time,psi)...)
            end
        end
        tout, ψ = timeevolution.schroedinger_dynamic(times, psi, get_hamiltonian,fout=feval)
        return (ψ[end][1], [[ψ[i][j] for i in 1:length(times)] for j in 2:length(ψ[1])]...)
    end
end
#,alg=RK4(),dt=(times[2]-times[1])/2,adaptive=false
    