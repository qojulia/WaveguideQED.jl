using Test
using WaveguideQED
using QuantumOpticsBase
using LinearAlgebra

ξ(t, σ, t0) = sqrt(2 / σ) * (log(2) / pi)^(1 / 4) * exp(-2 * log(2) * (t - t0)^2 / σ^2)
ξ2(t1, t2, σ, t0) = ξ(t1, σ, t0) * ξ(t2, σ, t0)
ex(args...; kw...) = waveguide_evolution(args...; tol=1e-11, kw...)   # tight tolerance

# Reference that shares no code with the Krylov propagation: exp(-iHτ)ψ summed as a Taylor series on
# `substeps` equal pieces of each bin (short enough that the series converges without cancellation).
# A time-dependent H is evaluated at the midpoint of each piece. Same output conventions as
# waveguide_evolution.
function taylor_evolution(times, ψ0, H; fout=nothing, substeps=10)
    ops = get_waveguide_operators(H)
    dt = get_dt(H.basis_l)
    times_sim = 0:dt:times[end]+dt
    τ = dt / substeps
    ψ, term, tmp = copy(ψ0), copy(ψ0), copy(ψ0)
    fvals = []
    states = [copy(ψ0)]
    for k in 1:length(times_sim)-1
        fout isa Function && push!(fvals, fout(times_sim[k], ψ))
        for j in 1:substeps
            H isa TimeDependentSum && set_time!(H, times_sim[k] + (j - 1 / 2) * τ)
            set_waveguidetimeindex!(ops, k)   # after fout, which may move shared operators
            term.data .= ψ.data
            for n in 1:100
                mul!(tmp, H, term, true, false)
                tmp.data .*= -im * τ / n
                term, tmp = tmp, term
                ψ.data .+= term.data
                norm(term.data) <= 1e-17 * norm(ψ.data) && break
            end
        end
        push!(states, copy(ψ))
    end
    fout === nothing && return ψ
    fout == 1 && return states
    (ψ, [getindex.(fvals, j) for j in eachindex(first(fvals))]...)
end

@testset "output conventions and exactness per bin" begin
    times = 0:0.1:2
    dt = 0.1
    bw = WaveguideBasis(1, 1, times)
    bc = FockBasis(1)
    H = im * sqrt(1 / dt) * (create(bc) ⊗ destroy(bw) - destroy(bc) ⊗ create(bw))
    ψ0 = fockstate(bc, 1) ⊗ zerophoton(bw)
    n = (create(bc) * destroy(bc)) ⊗ identityoperator(bw)
    fout(t, ψ) = (t, real(expect(n, ψ)))
    ψe, te, ne = waveguide_evolution(times, ψ0, H; fout=fout)
    ψr, tr, nr = taylor_evolution(times, ψ0, H; fout=fout)
    @test te == tr == collect(times)
    # every bin is an exact rotation between the cavity and the empty bin: n(k dt) = cos²ᵏ(√(γ dt))
    @test ne ≈ cos(sqrt(dt)) .^ (2 .* (0:length(times)-1)) atol = 1e-12
    @test nr ≈ ne atol = 1e-9
    @test ψe.data ≈ ψr.data atol = 1e-9
    @test waveguide_evolution(times, ψ0, H).data ≈ ψe.data
    states_e = waveguide_evolution(times, ψ0, H; fout=1)
    states_r = taylor_evolution(times, ψ0, H; fout=1)
    @test length(states_e) == length(states_r) == length(times) + 1
    @test all(isapprox(a.data, b.data; atol=1e-9) for (a, b) in zip(states_e, states_r))
    @test_logs (:warn, r"ignores the keyword arguments") waveguide_evolution(times, ψ0, H; abstol=1e-8)
    if Base.get_extension(WaveguideQED, :WaveguideQEDQuantumOpticsExt) === nothing
        err = try
            waveguide_montecarlo(times, ψ0, H, [n])
        catch e
            e
        end
        @test err isa MethodError
        @test occursin("using QuantumOptics", sprint(showerror, err))
    end
end

@testset "two photons on a cavity" begin
    times = 0:0.1:6
    dt = 0.1
    bw = WaveguideBasis(2, 1, times)
    bc = FockBasis(2)
    H = im * sqrt(1 / dt) * (create(bc) ⊗ destroy(bw) - destroy(bc) ⊗ create(bw)) +
        0.5 * (create(bc) * create(bc) * destroy(bc) * destroy(bc)) ⊗ identityoperator(bw)
    ψ0 = fockstate(bc, 0) ⊗ (twophoton(bw, ξ2, 1, 3) / sqrt(2))
    ψe = ex(times, ψ0, H)
    @test norm(ψe) ≈ 1 atol = 1e-8
    @test waveguide_evolution(times, ψ0, H).data ≈ ψe.data atol = 1e-6   # default tolerance
    @test ψe.data ≈ taylor_evolution(times, ψ0, H).data atol = 1e-8
    # a bin needs 6 Krylov vectors here, so maxiter=5 forces the steps to be split, which must not change
    # the result (smaller maxiter also works but is slow: the error estimate then shrinks only like τ^(maxiter-1))
    @test ex(times, ψ0, H; maxiter=5).data ≈ ψe.data atol = 1e-8
end

@testset "beamsplitter and swap between waveguides" begin
    times = 0:0.1:6
    dt = 0.1
    bw = WaveguideBasis(2, 2, times)
    w1, wd1, w2, wd2 = destroy(bw, 1), create(bw, 1), destroy(bw, 2), create(bw, 2)
    ψ0 = twophoton(bw, [1, 2], ξ2, 1, 3)
    for V in (pi / 4, pi / 2)
        H = im * V / dt * (wd2 * w1 - wd1 * w2)
        ψe = ex(times, ψ0, H)
        @test ψe.data ≈ taylor_evolution(times, ψ0, H).data atol = 1e-8
        V == pi / 4 && @test norm(TwoPhotonView(ψe, 1))^2 ≈ 0.5 atol = 1e-6   # Hong-Ou-Mandel bunching
    end
end

@testset "delayed feedback with operators shared by fout" begin
    tau = 1
    dt = 0.1
    times = 0:dt:2 * tau + dt
    bw = WaveguideBasis(1, 1, times)
    be = FockBasis(1)
    w, wd = destroy(bw), create(bw)
    w_tau, wd_tau = destroy(bw; delay=tau / dt + 1), create(bw; delay=tau / dt + 1)
    sd, s = create(be), destroy(be)
    Ie, Iw = identityoperator(be), identityoperator(bw)
    H = im * sqrt(1 / dt) * (sd ⊗ Ie ⊗ w - s ⊗ Ie ⊗ wd) + im * sqrt(1 / dt) * (Ie ⊗ sd ⊗ w_tau - Ie ⊗ s ⊗ wd_tau)
    nw = Ie ⊗ Ie ⊗ (wd * w)            # shares w, wd with H; expect_waveguide moves their time index
    ne_a = (sd * s) ⊗ Ie ⊗ Iw
    f(t, ψ) = (expect(ne_a, ψ), expect_waveguide(nw, ψ))
    ψ0 = fockstate(be, 1) ⊗ fockstate(be, 0) ⊗ zerophoton(bw)
    ψe, ae, we = ex(0:dt:4 * tau, ψ0, H; fout=f)
    ψr, ar, wr = taylor_evolution(0:dt:4 * tau, ψ0, H; fout=f)
    @test ψe.data ≈ ψr.data atol = 1e-8
    @test ae ≈ ar atol = 1e-8
    @test we ≈ wr atol = 1e-8
end

@testset "time-dependent Hamiltonian with piecewise-constant controls" begin
    times = 0:0.1:4
    dt = 0.1
    bw = WaveguideBasis(1, 1, times)
    bc = FockBasis(1)
    κ(t) = t < 1.23 ? 0.2 : (t < 2.71 ? 3.0 : 1.0)      # switching times inside bins
    δ(t) = t < 1.23 ? 0.0 : 2.0
    H = TimeDependentSum([t -> sqrt(κ(t) / dt), δ],
                         [im * (create(bc) ⊗ destroy(bw) - destroy(bc) ⊗ create(bw)),
                          (create(bc) * destroy(bc)) ⊗ identityoperator(bw)])
    ψ0 = fockstate(bc, 1) ⊗ zerophoton(bw)
    ψe = ex(times, ψ0, H; tstops=[1.23, 2.71])
    ψr = taylor_evolution(times, ψ0, H)   # the switching times fall on boundaries of the pieces
    @test ψe.data ≈ ψr.data atol = 1e-8
    # without tstops the switch inside a bin is only resolved to the midpoint rule
    @test !isapprox(ex(times, ψ0, H).data, ψr.data; atol=1e-6)
end

@testset "non-Hermitian Hamiltonian (Arnoldi)" begin
    times = 0:0.1:3
    dt = 0.1
    bw = WaveguideBasis(1, 1, times)
    bc = FockBasis(1)
    H = im * sqrt(1 / dt) * (create(bc) ⊗ destroy(bw) - destroy(bc) ⊗ create(bw)) -
        0.3im * (create(bc) * destroy(bc)) ⊗ identityoperator(bw)
    ψ0 = fockstate(bc, 1) ⊗ zerophoton(bw)
    ψe = ex(times, ψ0, H)
    @test ψe.data ≈ taylor_evolution(times, ψ0, H).data atol = 1e-8
    @test norm(ψe) < 1
end
