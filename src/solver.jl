"""
    waveguide_evolution(times, psi0, H; fout=nothing, tol=1e-8, maxiter=30, tstops=(), substeps=1, hermitian=nothing)

Evolve `psi0` under the time-binned Hamiltonian `H`. The waveguide operators in `H` address one time
bin at a time, so `H` is constant within each bin of width `dt` and changes discontinuously at the
bin boundaries. Each bin is therefore propagated by applying ``\\exp(-iH\\Delta t)`` directly with an
error-controlled Krylov method (Lanczos for Hermitian `H`, Arnoldi otherwise). The evolution starts at
`t=0` and covers the bins `0:dt:times[end]+dt`.

# Arguments
* `times`: Points of time for which output should be displayed (`fout` is evaluated at `0:dt:times[end]`).
* `psi0`: Initial state vector can only be a ket.
* `H`: Operator containing a [`WaveguideOperator`](@ref) either through a LazySum or LazyTensor. Time-dependent
  operators (e.g. a `TimeDependentSum`) are supported; their clock is set with `set_time!`.
* `fout=nothing`: If given, this function `fout(t, psi)` is called every time step. Example: `fout(t,psi) = expect(A,psi)` will return the epectation value of A at everytimestep.
   If `fout =1` the state psi is returned for all timesteps in a vector.
   ATTENTION: The state `psi` is neither normalized nor permanent! It is still in use by the solver and therefore must not be changed.
* `tol=1e-8`: Error per bin relative to the norm of the state.
* `maxiter=30`: Maximal Krylov dimension; steps that need more are split automatically.
* `tstops=()`: Times inside bins where a time-dependent `H` jumps, e.g. switching times of piecewise-constant controls.
* `substeps=1`: Split each bin into equal pieces; the coefficients of a time-dependent `H` are evaluated at the midpoint of each piece.
* `hermitian=nothing`: Whether `H` is Hermitian (`nothing` = detect automatically).

# Returns
* if `fout=nothing` the output of the solver will be the state `ψ` at the last timestep.
* if `fout` is given a tuple with the state `ψ` at the last timestep and the output of `fout` is given. If `fout` returns a tuple the tuple will be flattened.
* if `fout = 1` `ψ` at all timesteps is returned.

# Examples

* `fout(t,psi) = (expect(A,psi),expect(B,psi))` will result in  a tuple (ψ, ⟨A(t)⟩,⟨B(t)⟩), where `⟨A(t)⟩` is a vector with the expectation value of `A` as a function of time.

"""
function waveguide_evolution(times, psi, H; fout=nothing, tol=1e-8, maxiter::Int=30, tstops=(),
                             substeps::Int=1, hermitian=nothing, kwargs...)
    isempty(kwargs) || @warn "waveguide_evolution ignores the keyword arguments $(Tuple(keys(kwargs))). Each time bin is propagated exactly with a Krylov method (no ODE solver); its accuracy is set with `tol`." maxlog=1
    ops = get_waveguide_operators(H)
    dt = get_dt(H.basis_l)
    tend = times[end]+dt
    times_sim = 0:dt:tend
    nbins = length(times_sim) - 1
    _check_initial_norm(psi)

    time_dependent = !QuantumOpticsBase.is_const(H)
    cuts = sort!([Float64(t) for t in tstops if 0 < t < tend])
    function set_bin!(k, t)
        time_dependent && set_time!(H, t)
        set_waveguidetimeindex!(ops, k)
    end
    set_bin!(1, dt/2)
    herm = hermitian === nothing ? _probably_hermitian(H, psi) : hermitian
    cache = KrylovCache(psi, maxiter, herm)
    # tolerances below the floating point resolution of the state can never be met
    tol_eff = max(tol, 100*eps(real(eltype(psi.data))))

    ψ = copy(psi)
    out = copy(psi)
    fvals = fout === nothing || fout == 1 ? nothing : Any[]
    states = fout == 1 ? [copy(psi)] : nothing
    edges = Float64[]
    for k in 1:nbins
        t0, t1 = times_sim[k], times_sim[k+1]
        fvals === nothing || push!(fvals, (fout(t0, ψ)...,))
        empty!(edges)
        push!(edges, t0)
        for t in cuts
            t0 + 1e-9*dt < t < t1 - 1e-9*dt && push!(edges, t)
        end
        push!(edges, t1)
        for s in 1:length(edges)-1, j in 1:substeps
            ta = edges[s] + (j-1)*(edges[s+1]-edges[s])/substeps
            tb = edges[s] + j*(edges[s+1]-edges[s])/substeps
            set_bin!(k, (ta + tb)/2)
            _expv_step!(out, H, ψ, tb - ta, cache, tol_eff)
            ψ, out = out, ψ
        end
        states === nothing || push!(states, copy(ψ))
    end
    states === nothing || return states
    _evolution_output(ψ, fvals)
end

function _check_initial_norm(psi)
    isapprox(norm(psi),1,rtol=10^(-6)) || @warn "Initial waveguidestate is not normalized. Consider passing norm=true to the state generation function."
end

# Output as returned by waveguide_evolution: `fvals[k]` holds the (splatted) output of fout at the k'th time.
function _evolution_output(ψ, fvals)
    fvals === nothing && return ψ
    (ψ, [[f[j] for f in fvals] for j in eachindex(first(fvals))]...)
end

# Workspace for Krylov propagation ψ ↦ exp(-iHτ)ψ. Krylov vectors are allocated on first use, so memory
# grows only to the subspace dimension that is actually needed.
mutable struct KrylovCache{K,T}
    V::Vector{K}              # Krylov basis
    w::K                      # H*v_j
    tmp::Vector{K}            # intermediate states when a step has to be split
    maxiter::Int
    hermitian::Bool
    α::Vector{Float64}        # Lanczos: diagonal of the tridiagonal projection
    β::Vector{Float64}        # Lanczos: off-diagonal
    Hm::Matrix{ComplexF64}    # Arnoldi: Hessenberg projection
    nmul::Int                 # number of H*ψ products (diagnostics)
end
function KrylovCache(psi::Ket, maxiter::Int, hermitian::Bool)
    maxiter >= 1 || throw(ArgumentError("maxiter must be positive"))
    KrylovCache{typeof(psi),eltype(psi.data)}(typeof(psi)[], copy(psi), typeof(psi)[], maxiter, hermitian,
        zeros(maxiter), zeros(maxiter), zeros(ComplexF64, maxiter+1, maxiter), 0)
end
function _krylov_vector!(c::KrylovCache, j, template)
    while length(c.V) < j
        push!(c.V, Ket(template.basis, similar(template.data)))
    end
    c.V[j]
end

# ψ ↦ exp(-iHτ)ψ; steps where the Krylov iteration does not converge within `maxiter` are halved.
function _expv_step!(out, H, ψ, τ, c::KrylovCache, tol, depth=0)
    _krylov_expv!(out, H, ψ, τ, c, tol) && return out
    depth < 40 || error("Krylov propagation did not converge; increase `maxiter` or `tol`.")
    length(c.tmp) <= depth && push!(c.tmp, copy(ψ))
    mid = c.tmp[depth+1]
    _expv_step!(mid, H, ψ, τ/2, c, tol, depth+1)
    _expv_step!(out, H, mid, τ/2, c, tol, depth+1)
end

# exp(-iτT)e₁ for the symmetric tridiagonal T = tridiag(β, α, β)
function _expT_e1(α, β, j, τ)
    F = eigen(SymTridiagonal(α[1:j], β[1:j-1]))
    F.vectors * (cis.(-τ .* F.values) .* F.vectors[1, :])
end

# Krylov approximation of exp(-iHτ)ψ, written to `out`. The iteration stops when Saad's a-posteriori
# error estimate `β_j |eⱼᵀexp(-iτT)e₁|` drops below `tol` (relative to `norm(ψ)`). Returns `false` if
# that did not happen within `c.maxiter` iterations.
function _krylov_expv!(out::Ket, H, ψ::Ket, τ::Real, c::KrylovCache{K,T}, tol) where {K,T}
    R = real(T)
    β0 = norm(ψ.data)
    if iszero(β0)
        fill!(out.data, zero(T))
        return true
    end
    v = _krylov_vector!(c, 1, ψ)
    v.data .= ψ.data .* R(inv(β0))
    w = c.w
    for j in 1:c.maxiter
        vj = c.V[j]
        mul!(w, H, vj, true, false)
        c.nmul += 1
        if c.hermitian
            α = real(dot(vj.data, w.data))
            c.α[j] = α
            βj = _lanczos_orthogonalize!(w.data, vj.data, j > 1 ? c.V[j-1].data : nothing, R(α), R(j > 1 ? c.β[j-1] : 0))
            c.β[j] = βj
            y = _expT_e1(c.α, c.β, j, τ)
        else
            for i in 1:j
                h = dot(c.V[i].data, w.data)
                c.Hm[i, j] = h
                axpy!(T(-h), c.V[i].data, w.data)
            end
            βj = norm(w.data)
            c.Hm[j+1, j] = βj
            y = exp(-im * τ * c.Hm[1:j, 1:j])[:, 1]
        end
        if βj * abs(y[j]) <= tol
            _krylov_combine!(out.data, c.V, [T(β0 * y[i]) for i in 1:j])
            return true
        end
        j == c.maxiter && return false
        vn = _krylov_vector!(c, j+1, ψ)
        vn.data .= w.data .* R(inv(βj))
    end
    return false
end

# w ← w - α v - β u (u = previous Lanczos vector or nothing); returns norm(w). Fused into a single
# pass for CPU vectors, since for cheap Hamiltonians the Krylov vector operations dominate.
function _lanczos_orthogonalize!(w::Vector{T}, v::Vector{T}, u, α::Real, β::Real) where {T}
    s = zero(real(T))
    if u === nothing
        @inbounds @simd for i in eachindex(w, v)
            x = w[i] - α * v[i]
            w[i] = x
            s += abs2(x)
        end
    else
        @inbounds @simd for i in eachindex(w, v, u)
            x = w[i] - α * v[i] - β * u[i]
            w[i] = x
            s += abs2(x)
        end
    end
    sqrt(s)
end
function _lanczos_orthogonalize!(w, v, u, α, β)
    axpy!(-α, v, w)
    u === nothing || axpy!(-β, u, w)
    norm(w)
end

# out = Σᵢ coef[i] V[i] (single pass for CPU vectors)
function _krylov_combine!(out::Vector{T}, V, coef) where {T}
    Vd = [V[i].data for i in eachindex(coef)]
    @inbounds for k in eachindex(out)
        s = zero(T)
        for i in eachindex(coef)
            s += coef[i] * Vd[i][k]
        end
        out[k] = s
    end
    out
end
function _krylov_combine!(out, V, coef)
    fill!(out, zero(eltype(out)))
    for i in eachindex(coef)
        axpy!(coef[i], V[i].data, out)
    end
    out
end

# Test ⟨x|Hy⟩ = conj(⟨y|Hx⟩) at the current time/bin of H with two generic (deterministic, so the
# global RNG is left alone) vectors.
function _probably_hermitian(H, psi::Ket)
    T = eltype(psi.data)
    n = length(psi.data)
    x, y, Hx, Hy = (copy(psi) for _ in 1:4)
    copyto!(x.data, T[cis(sqrt(2) * i^1.3) * (1 + i / n) for i in 1:n])
    copyto!(y.data, T[cis(sqrt(3) * i^1.1) * (2 - i / n) for i in 1:n])
    mul!(Hx, H, x, true, false)
    mul!(Hy, H, y, true, false)
    a = dot(x.data, Hy.data)
    b = conj(dot(y.data, Hx.data))
    abs(a - b) <= sqrt(eps(real(T))) * (norm(Hx.data) * norm(y.data) + norm(Hy.data) * norm(x.data))
end

"""
    waveguide_montecarlo(times,psi,H,J;fout=nothing)

See documentation for [`waveguide_evolution`](@ref) on how to define `fout`. J should be a list of collapse operators following documentation of [`timeevolution.mcwf_dynamic`](https://docs.qojulia.org/api/#QuantumOptics.timeevolution.mcwf_dynamic).

Requires [`QuantumOptics.jl`](https://qojulia.org/), which provides the solver: the method is defined once `using QuantumOptics` has been run.
"""
function waveguide_montecarlo end

"""
    fast_unitary(times_eval,psi,H;order=2,fout=nothing)

See documentation for [`waveguide_evolution`](@ref) on how to define `fout`. J should be a list of collapse operators following documentation of [`timeevolution.mcwf_dynamic`](https://docs.qojulia.org/api/#QuantumOptics.timeevolution.mcwf_dynamic). 
"""
function fast_unitary(times_eval,psi,H;order=2,fout=nothing)
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
