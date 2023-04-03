abstract type LazyKet{B} end

"""
    LazyTensorKet(kets)

Lazy tensor product between kets. Used in functions for beamsplitter and subsequent detection.
"""
mutable struct LazyTensorKet{B} <: LazyKet{B}
    basis::B
    kets
    factor
end
function LazyTensorKet(kets::Ket...)
    btotal = basis(kets[1]) 
    for ψ in kets[2:end]
        btotal = btotal ⊗ basis(ψ)
    end
    B = typeof(btotal)
    LazyTensorKet{B}(btotal,collect(kets),1)
end
Base.copy(x::LazyTensorKet) = LazyTensorKet(deepcopy(x.kets)...)
function Base.:*(a::Number,x::LazyTensorKet)
    out = copy(x)
    out.factor *= a
    out
end
Base.:/(x::LazyTensorKet,a::Number) = *(1/a,x)
dagger(x::LazyTensorKet) = LazyTensorBra(x.basis,dagger.(x.kets),x.factor)

"""
    LazyTensorBra(bras)

Lazy tensor product between bras. Used in functions for beamsplitter and subsequent detection.
"""
mutable struct LazyTensorBra{B} <: LazyKet{B}
    basis::B
    bras
    factor
end
function LazyTensorBra(bras::Bra...)
    btotal = basis(bras[1]) 
    for ψ in bras[2:end]
        btotal = btotal ⊗ basis(ψ)
    end
    B = typeof(btotal)
    LazyTensorBra{B}(btotal,collect(bras),1)
end
Base.copy(x::LazyTensorBra) = LazyTensorBra(deepcopy(x.bras))
function Base.:*(a::Number,x::LazyTensorBra)
    out = copy(x)
    out.factor *= a
    out
end
Base.:/(x::LazyTensorBra,a::Number) = *(1/a,x)
dagger(x::LazyTensorBra) = LazyTensorKet(x.basis,dagger.(x.bras),x.factor)

"""
    LazySumKet(kets...)

Lazy sum of LazyTensorKets that is used to express entanglement between subsystems in LazyTensorKets. 
"""
mutable struct LazySumKet{B,F} <: LazyKet{B}
    basis::B
    kets
    factors::F
end
function LazySumKet(kets,factors) 
    for k in kets[2:end]
        @assert kets[1].basis == k.basis
    end
    @assert length(kets) == length(factors)
    LazySumKet(kets[1].basis,kets,factors)
end
LazySumKet(kets::LazyTensorKet...) = LazySumKet([kets...],ones(ComplexF64, length(kets)))
Base.copy(x::LazySumKet) = LazySumKet(x.basis,copy(x.kets),x.factors)
function Base.:*(a::Number,x::LazySumKet)
    out = copy(x)
    out.factors *= a
    out
end
Base.:/(x::LazySumKet,a::Number) = *(1/a,x)
Base.:+(a::LazySumKet{B,F}, b::LazySumKet{B,F}) where {B,F} = LazySumKet([a.kets..., b.kets...],[a.factors; b.factors])
Base.:-(a::LazySumKet{B,F}, b::LazySumKet{B,F}) where {B,F} = LazySumKet([a.kets..., b.kets...],[a.factors; -b.factors])
function Base.:+(a::LazyTensorKet,b::LazyTensorKet)
    @assert a.basis==b.basis
    LazySumKet(a)+LazySumKet(b)
end
function Base.:-(a::LazyTensorKet,b::LazyTensorKet)
    @assert a.basis==b.basis
    LazySumKet(a) - LazySumKet(b)
end
function Base.:+(a::LazySumKet,b::LazyTensorKet)
    @assert a.basis==b.basis
    a + LazySumKet(b)
end
function Base.:+(a::LazyTensorKet,b::LazySumKet)
    @assert a.basis==b.basis
    LazySumKet(a) + b
end
function Base.:-(a::LazySumKet,b::LazyTensorKet)
    @assert a.basis==b.basis
    a - LazySumKet(b)
end
function Base.:-(a::LazyTensorKet,b::LazySumKet)
    @assert a.basis==b.basis
    LazySumKet(a) - b
end
function Base.:*(a::LazyTensorKet,b::LazyTensorKet)
    @assert a.basis==b.basis
    a.factor*b.factor*(dagger.(a.kets).*b.kets)
end
function Base.:*(a::LazyTensorBra,b::LazyTensorKet)
    @assert a.basis==b.basis
    a.factor*b.factor*(a.bras.*b.kets)
end
get_tmp(x::LazyTensorKet) = copy(x)
get_tmp(x::LazySumKet) = copy(x.kets[1])

"""
    Detector(wa,wb)

Detector operation defined by giving waveguide annihilation operator `wa` and `wb` from two subsystems. `wa` acts on the first subsystem of a [`LazyTensorKet`](@ref) or [`LazySumKet`](@ref) and `wb` on the second.
"""
mutable struct Detector{B}
    basis::B
    operators
    function Detector(wa,wb)
        btotal = wa.basis_l ⊗ wb.basis_l
        B = typeof(btotal)
        new{B}(btotal,[wa,wb])
    end
end
set_waveguidetimeindex!(detector::Detector,index) = set_waveguidetimeindex!(detector.operators,index)

function mul!(result::LazyTensorKet{B},detector::Detector{B},input::LazyTensorKet{B},alpha,beta) where {B}
    for i in eachindex(detector.operators)
        mul!(result.kets[i],detector.operators[i],input.kets[i],alpha,beta)
    end
end
function mul!(result::LazySumKet{B,F},detector::Detector{B},input::LazySumKet{B,F},alpha,beta) where {B,F}
    for i in eachindex(result.kets)
        mul!(result.kets[i],detector,input.kets[i],alpha,beta)
    end
end
function get_waveguide_operators(x::Detector)
    get_waveguide_operators(x.operators)
end

"""
    detect_single_click(ψ,detector::Detector,projection)
    detect_single_click(ψ,detector::Detector)

Calculate probability of observing `projection` after beamsplitter operation and subsequent detection event defined in `detector` on the state `ψ`.

# Arguments
* `ψ` can be either [`LazyTensorKet`](@ref) or [`LazySumKet`](@ref) and is the state on which the beamsplitter and detection is applied
* `detector` defines the beamsplitter and subsequent detection operation. See [`Detector`](@ref) for more details on how to define.
* `projection` if given is a [`LazyTensorKet`](@ref) or [`LazySumKet`](@ref) which projects onto the state after the measurement.  If no projection is given, instead the total probability of having the detector click is given by applying all possible combinations of projections using [`get_all_projectors`](@ref).

# Returns
* If `projection` is given returns probability of having detector click and being in state defined by `projection`
* If `projection` is not given returns the total probability of having a the detector click (only a single click, for double clicks use [`detect_double_click`](@ref)) by applying all possibile projections with zerophotons in the waveguide using [`get_all_projectors`](@ref).

See [Beamsplitter](https://mabuni1998.github.io/WaveguideQED/dev/detection_example/) for an example on how to use.
"""
function detect_single_click(ψ,detector::Detector,projection)
    waveguide_operators=get_waveguide_operators(detector)
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    p_time = zeros(ComplexF64,timesteps)
    tmp  = get_tmp(ψ)
    detect_single_click!(p_time,tmp,ψ,detector,projection)
    p_time,real(transpose(conj.(p_time))*p_time)
end
function detect_single_click(ψ,detector::Detector)
    p1s = get_all_projectors(detector.operators[1].basis_l)
    p2s = get_all_projectors(detector.operators[1].basis_l)
    total_detect = 0
    waveguide_operators=get_waveguide_operators(detector)
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    p_time = zeros(ComplexF64,timesteps)
    tmp  = get_tmp(ψ)
    for p1 in p1s
        for p2 in p2s
            p_time.=0
            detect_single_click!(p_time,tmp,ψ,detector,LazyTensorKet(p1,p2))
            total_detect += real(transpose(conj.(p_time))*p_time)
        end
    end
    total_detect
end
function detect_single_click!(p_time,detected_state::LazyTensorKet{B},ψ::LazyTensorKet{B},detector::Detector{B},projection::LazyTensorKet{B};timeindex = 1,factor=1) where {B}
    waveguide_operators=get_waveguide_operators(detector)
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    projector = dagger(projection)
    proj_undetected_a,proj_undetected_b = projector*ψ
    if iszero(proj_undetected_a) & iszero(proj_undetected_b)
        return p_time
    end 
    for i in timeindex:timesteps
        set_waveguidetimeindex!(detector,i)
        mul!(detected_state,detector,ψ,1,0)
        proj_detected_a,proj_detected_b = projector*detected_state
        p_time[i] += factor*proj_undetected_a*proj_detected_b + factor*proj_detected_a*proj_undetected_b
    end
    p_time
end
function detect_single_click!(p_time,detected_state::LazyTensorKet{B},ψ::LazySumKet{B,F},detector::Detector{B},projection::LazyTensorKet{B};timeindex=1,factor=1) where {B,F}
    for (i,psi) in enumerate(ψ.kets)
        detect_single_click!(p_time,detected_state,psi,detector,projection,timeindex=timeindex,factor = factor*ψ.factors[i])
    end
    p_time
end
function detect_single_click!(p_time,detected_state::LazyTensorKet{B},ψ::LazySumKet{B,F},detector::Detector{B},projection::LazySumKet{B,F};timeindex=1,factor=1) where {B,F}
    for (i,psi) in enumerate(ψ.kets)
        detect_single_click!(p_time,detected_state,psi,detector,projection,timeindex=timeindex,factor=factor*ψ.factors[i])
    end
    p_time
end
function detect_single_click!(p_time,detected_state::LazyTensorKet{B},ψ::LazyTensorKet{B},detector::Detector{B},projection::LazySumKet{B,F};timeindex = 1,factor=1) where {B,F}
    waveguide_operators=get_waveguide_operators(detector)
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    @assert timesteps == length(p_time)
    for i in timeindex:timesteps
        set_waveguidetimeindex!(detector,i)
        mul!(detected_state,detector,ψ,1,0)
        for (j,p) in enumerate(projection.kets)
            proj_undetected_a,proj_undetected_b = p*ψ
            proj_detected_a,proj_detected_b = p*detected_state
            p_time[i] += projection.factors[j]*factor*(proj_undetected_a*proj_detected_b + proj_detected_a*proj_undetected_b)
        end
    end
    p_time
end

Base.:*(projection::LazyKet{B},detector::Detector{B},ψ::LazyKet{B}) where {B} = detect_single_click(ψ,detector,projection)
Base.:*(detector::Detector{B},ψ::LazyKet{B}) where {B} = detect_single_click(ψ,detector)


"""
    get_all_projectors(b)

Returns all combinations of possible states with zerophotons in the waveguide. 
"""
function get_all_projectors(b::CompositeBasis)
    projectors = Array{Any}(undef,length(b.bases))
    for (i,basis) in enumerate(b.bases)
        if isa(basis,WaveguideBasis)
            projectors[i] = [zerophoton(basis)]
        else
            projectors[i] = [basisstate(basis,j) for j in 1:length(basis)]
        end
    end
    out = projectors[1]
    for p in projectors[2:end]
       out = mul_internal(out,p) 
    end
    out
end
function get_all_projectors(b::WaveguideBasis)
    [zerophoton(b)]
end
function mul_internal(p1,p2)
    out = Array{Ket}(undef,length(p1)*length(p2))
    for i in eachindex(p1)
        for j in eachindex(p2)
            out[(i-1)*length(p2)+j] = p1[i] ⊗ p2[j]
        end 
    end
    out
end



function _two_time_total_probability(p_time,timesteps)
    out = 0
    @simd for i in 1:timesteps
        @simd for j in i:timesteps
            @inbounds out += (i==j ? 1/2 : 1) * real(p_time[i,j] * conj(p_time[i,j]))
        end
    end
    out
end

"""
    detect_double_click(ψ,detector1,detector2,projection)
    detect_double_click(ψ,detector1,detector2)

Calculate probability of observing `projection` after beamsplitter operation and two subsequent detection events defined in `detector1` and `detector2` on the state `ψ`.

# Arguments
* `ψ` can be either [`LazyTensorKet`](@ref) or [`LazySumKet`](@ref) and is the state on which the beamsplitter and detection is applied
* `detector1` defines the first beamsplitter and subsequent detection operation. See [`Detector`](@ref) for more details on how to define.
* `detector2` defines the second beamsplitter and subsequent detection operation. See [`Detector`](@ref) for more details on how to define.
* `projection` if given is a [`LazyTensorKet`](@ref) or [`LazySumKet`](@ref) which projects onto the state after the measurement.  If no projection is given, instead the total probability of having the detector click is given by applying all possible combinations of projections using [`get_all_projectors`](@ref).

# Returns
* If `projection` is given: Returns probability of having `detector1` and `detector2` click and being in state defined by `projection`
* If `projection` is not given: Returns the total probability of having `detector1` and `detector2` click by applying all possibile projections with zerophotons in the waveguide using [`get_all_projectors`](@ref).

See [Beamsplitter](https://mabuni1998.github.io/WaveguideQED/dev/detection_example/) for an example on how to use. 

"""
function detect_double_click(ψ,detector1::Detector{B},detector2::Detector{B},projection) where {B}
    waveguide_operators=[get_waveguide_operators(detector1)...,get_waveguide_operators(detector2)...]
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    p_time = zeros(ComplexF64,(timesteps,timesteps))
    a_measured  = get_tmp(ψ)
    tmp  = get_tmp(ψ)
    detect_double_click!(p_time,a_measured,tmp,ψ,detector1,detector2,projection)
    p_time,_two_time_total_probability(p_time,timesteps)
end
function detect_double_click(ψ,detector1::Detector{B},detector2::Detector{B}) where B
    p1s = get_all_projectors(detector1.operators[1].basis_l)
    p2s = get_all_projectors(detector1.operators[2].basis_l)
    total_detect = 0
    waveguide_operators=[get_waveguide_operators(detector1)...,get_waveguide_operators(detector2)...]
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    p_time = zeros(ComplexF64,(timesteps,timesteps))
    a_measured  = get_tmp(ψ)
    tmp  = get_tmp(ψ)
    for p1 in p1s
        for p2 in p2s
            p_time.=0
            detect_double_click!(p_time,a_measured,tmp,ψ,detector1,detector2,LazyTensorKet(p1,p2))
            total_detect += _two_time_total_probability(p_time,timesteps)
        end
    end
    total_detect
end
function detect_double_click!(p_time,measure1::LazyTensorKet{B},tmp::LazyTensorKet{B},ψ::LazyTensorKet{B},detector1::Detector{B},detector2::Detector{B},projection::LazyTensorKet{B};factor=1,timeindex=1) where {B}
    waveguide_operators=[get_waveguide_operators(detector1)...,get_waveguide_operators(detector2)...]
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    @assert timesteps == size(p_time)[1]
    @assert timesteps == size(p_time)[2]
    for i in timeindex:timesteps
        set_waveguidetimeindex!(detector1,i)
        mul!(measure1,detector1,ψ,1,0)
        measured_state = LazySumKet(LazyTensorKet(measure1.kets[1],ψ.kets[2]),LazyTensorKet(ψ.kets[1],measure1.kets[2]))
        detect_single_click!(view(p_time,i,:),tmp,measured_state,detector2,projection;timeindex=i,factor=factor)
    end
    p_time
end
function detect_double_click!(p_time,measure1::LazyTensorKet{B},tmp::LazyTensorKet,ψ::LazyTensorKet{B},detector1::Detector{B},detector2::Detector{B},projection::LazySumKet{B,F};factor=1,timeindex=1) where {B,F}
    waveguide_operators=[get_waveguide_operators(detector1)...,get_waveguide_operators(detector2)...]
    timesteps = min([w.basis_l.nsteps for w in waveguide_operators]...) 
    @assert timesteps == size(p_time)[1]
    @assert timesteps == size(p_time)[2]
    for i in timeindex:timesteps
        set_waveguidetimeindex!(detector1,i)
        mul!(measure1,detector1,ψ,1,0)
        measured_state = LazySumKet(LazyTensorKet(measure1.kets[1],ψ.kets[2]),LazyTensorKet(ψ.kets[1],measure1.kets[2]))
        detect_single_click!(view(p_time,i,:),tmp,measured_state,detector2,projection;timeindex=i,factor=factor)
    end
    p_time
end
function detect_double_click!(p_time,a_measured::LazyTensorKet{B},tmp::LazyTensorKet{B},ψ::LazySumKet{B,F},detector1::Detector{B},detector2::Detector{B},projection::LazyTensorKet{B};timeindex=1,factor=1) where {B,F}
    for (i,psi) in enumerate(ψ.kets)
        detect_double_click!(p_time,a_measured,tmp,psi,detector1,detector2,projection,timeindex=timeindex,factor = factor*ψ.factors[i])
    end
    p_time
end
function detect_double_click!(p_time,a_measured::LazyTensorKet{B},tmp::LazyTensorKet{B},ψ::LazySumKet{B,F},detector1::Detector{B},detector2::Detector{B},projection::LazySumKet{B,F};timeindex=1) where {B,F}
    for (i,psi) in enumerate(ψ.kets)
        detect_double_click!(p_time,a_measured,tmp,psi,detector1,detector2,projection,timeindex=timeindex,factor=ψ.factors[i])
    end
    p_time
end

Base.:*(projection::LazyKet{B},detector2::Detector{B},detector1::Detector{B},ψ::LazyKet{B}) where {B} = detect_double_click(ψ,detector1,detector2,projection)
Base.:*(detector2::Detector{B},detector1::Detector{B},ψ::LazyKet{B}) where {B} = detect_double_click(ψ,detector1,detector2)