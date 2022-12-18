function get_waveguide_basis(basis::CompositeBasis)
    for b in basis.bases
        if isa(b,WaveguideBasis)
            return b
        end
    end
    error("No waveguide operator used. Use timeevolution.schroedinger from QuantumOptics.jl instead")
end

function waveguide_evolution(times,psi,H;fout=nothing)
    bw = get_waveguide_basis(H.basis_l)
    dt = times[2] - times[1]
    tend = times[end]
    function get_hamiltonian(time,psi)
        bw.timeindex=floor(Int,time/dt)+1
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
        tout, ψ = timeevolution.schroedinger_dynamic(times, psi, get_hamiltonian,fout=fout)
        return ψ
    end
end