module WaveguideQEDQuantumOpticsExt

using WaveguideQED
import WaveguideQED: waveguide_montecarlo, _check_initial_norm, _waveguide_timeindex
using QuantumOptics

function waveguide_montecarlo(times,psi,H,J;fout=nothing,kwargs...)
    ops = get_waveguide_operators(H)
    dt = times[2] - times[1]
    tend = times[end]
    Jdagger = dagger.(J)

    _check_initial_norm(psi)

    function get_hamiltonian(time,psi)
        set_waveguidetimeindex!(ops,_waveguide_timeindex(time, dt, RoundUp))
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

end
