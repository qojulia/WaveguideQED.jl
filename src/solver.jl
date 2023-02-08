#Used to extract WaveguideBasis to change index in get_hamiltonian()
function get_waveguide_basis(basis::CompositeBasis)
    for b in basis.bases
        if isa(b,WaveguideBasis) || isa(b,OneTimeBasis)
            return b
        end
    end
    error("No waveguide operator used. Use timeevolution.schroedinger from QuantumOptics.jl instead")
end

#Call solver from QuantumOptics.jl by defining functions for extraction and updating Hamiltonian.
function waveguide_evolution(times,psi,H;fout=nothing)
    basis = get_waveguide_basis(psi.basis)
    dt = times[2] - times[1]
    tend = times[end]
    function get_hamiltonian(time,psi)
        basis.timeindex= round(Int,time/dt) +1
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
        tout, ψ = timeevolution.schroedinger_dynamic(times, psi, get_hamiltonian,fout=eval_last_element,dtmax=(times[2]-times[1])/2)
        return ψ[end]
    else
        function feval(time,psi)
            if time == tend
                return (psi,fout(time,psi)...)
            else
                return (0,fout(time,psi)...)
            end
        end
        tout, ψ = timeevolution.schroedinger_dynamic(times, psi, get_hamiltonian,fout=feval,dtmax=(times[2]-times[1])/2)
        return (ψ[end][1], [[ψ[i][j] for i in 1:length(times)] for j in 2:length(ψ[1])]...)
    end
end