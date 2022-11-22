export cwbasis,CWState,ZeroPhoton,OnePhoton,TwoPhoton

import Base:+,-,/,*,zero,copy
using QuantumOptics

  
#Basis of the CWState. Will contain more information in future
struct cwbasis
    times
    N_cutoff::Int
end


abstract type WaveGuideState end

"""
Structure containing the whole state of the CW

data contains list of state 
length(data) = basis.N_cutoff


"""
mutable struct CWState
    basis::cwbasis
    data::Vector{WaveGuideState}
    diff::Vector{WaveGuideState}
    tmp_data::Vector{ComplexF64}
    tmp_diff::Vector{ComplexF64}
    timeindex::Int
end

function mul!(state::CWState,a::ComplexF64)
    for i in 1:length(state.data)
        mul!(state.data[i],a,state.timeindex)
    end
end

function mul_diff!(state::CWState,a::ComplexF64)
    for i in 1:length(state.diff)
        mul!(state.diff[i],a,state.timeindex)
    end    
end

mul_diff!(state::CWState,a::Number) = mul_diff!(state,complex(float(a)))

mul!(a::ComplexF64,state::CWState) = mul!(state::CWState,a::ComplexF64)
mul!(a::Number,state::CWState) = mul!(state,complex(float(a)))
mul!(state::CWState,a::Number) = mul!(state,complex(float(a)))

function +(state1::CWState,state2::CWState)
    if state1.timeindex!=state2.timeindex
        print("Warning:Not same timeindex in +(CWState,CWState)")
    end
    for i in 1:length(state1.data)
        add!(state1.data[i],state2.data[i],state1.timeindex)
    end
end

function copy(state::CWState)
    CWState(state.basis,copy(state.data),copy(state.diff),copy(state.tmp_data),copy(state.tmp_diff),state.timeindex)
end

function copy(statedata::Vector{WaveGuideState})
    out = Array{WaveGuideState}(undef,length(statedata))
    for i in 1:length(statedata)
        out[i] = copy(statedata[i])
    end
    return out
end


function zero(b::cwbasis)
    if b.N_cutoff == 0
        CWState(b,[ZeroPhoton(complex(0.0))],[ZeroPhoton(complex(0.0))],zeros(ComplexF64,1),1,1)
    elseif b.N_cutoff == 1
        CWState(b,
        [ZeroPhoton(complex(0.0)),OnePhoton(complex(0.0),zeros(ComplexF64,length(b.times)))],
        [ZeroPhoton(complex(0.0)),OnePhoton(complex(0.0),zeros(ComplexF64,length(b.times)))],
        zeros(ComplexF64,2+length(b.times)),
        zeros(ComplexF64,2+length(b.times)),
        1)
    elseif b.N_cutoff == 2
        CWState(b,
        [ZeroPhoton(complex(0.0)),OnePhoton(complex(0.0),zeros(ComplexF64,length(b.times))),TwoPhoton(complex(0.0),zeros(ComplexF64,length(b.times)),zeros(ComplexF64,(length(b.times),length(b.times))))],
        [ZeroPhoton(complex(0.0)),OnePhoton(complex(0.0),zeros(ComplexF64,length(b.times))),TwoPhoton(complex(0.0),zeros(ComplexF64,length(b.times)),zeros(ComplexF64,(length(b.times),length(b.times))))],
        zeros(ComplexF64,3+3*length(b.times)),
        zeros(ComplexF64,3+3*length(b.times)),
        1) 
    end
end

function zero(state::CWState)
    zero(state.basis)
end


#Structures containing 0,1 or 2 photons in the cavity waveguide state
mutable struct ZeroPhoton <:WaveGuideState
    ξ0::ComplexF64
end

mutable struct OnePhoton <:WaveGuideState
    ξ10::ComplexF64
    ξ01::Vector{ComplexF64}
end

mutable struct TwoPhoton <:WaveGuideState
    ξ20::ComplexF64
    ξ11::Vector{ComplexF64}
    ξ02::Matrix{ComplexF64}
end

function copy(state::ZeroPhoton)
    return ZeroPhoton(state.ξ0)
end

function copy(state::OnePhoton)
    return OnePhoton(state.ξ10,copy(state.ξ01))
end

function copy(state::TwoPhoton)
    return TwoPhoton(state.ξ20,copy(state.ξ11),copy(state.ξ02))
end

function add!(state1::ZeroPhoton,state2::ZeroPhoton,timeindex::Int)
    state1.ξ0 = state1.ξ0 + state2.ξ0
end

function add!(state1::OnePhoton,state2::OnePhoton,timeindex::Int)
    state1.ξ10 = state1.ξ10 + state2.ξ10
    state1.ξ01 = state1.ξ01 + state2.ξ01
end
function add!(state1::TwoPhoton,state2::TwoPhoton,timeindex::Int)
    state1.ξ20 = state1.ξ20 + state2.ξ20
    state1.ξ11 = state1.ξ11 .+ state2.ξ11
    for j in 1:timeindex
        state1.ξ02[timeindex,j] = state1.ξ02[timeindex,j]+state2.ξ02[timeindex,j]
    end
    for j in timeindex+1:length(state1.ξ11)
        state1.ξ02[j,timeindex] = state1.ξ02[j,timeindex]+state2.ξ02[j,timeindex]
    end
end

function zero_diff!(state1::ZeroPhoton,timeindex::Int)
    state1.ξ0 = 0
end

function zero_diff!(state1::OnePhoton,timeindex::Int)
    state1.ξ10  = 0
    state1.ξ01 .= 0
end
function zero_diff!(state1::TwoPhoton,timeindex::Int)
    state1.ξ20 = 0
    state1.ξ11[1:end] .= 0
    for j in 1:timeindex
        state1.ξ02[timeindex,j] = 0
    end
    for j in timeindex+1:length(state1.ξ11)
        state1.ξ02[j,timeindex] = 0
    end
end




function set_equal!(state1::ZeroPhoton,state2::ZeroPhoton,timeindex::Int)
    state1.ξ0 = state2.ξ0
end

function set_equal!(state1::OnePhoton,state2::OnePhoton,timeindex::Int)
    state1.ξ10 = state2.ξ10
    state1.ξ01[1:end] = state2.ξ01[1:end]
end

function set_equal!(state1::TwoPhoton,state2::TwoPhoton,timeindex::Int)
    state1.ξ20 = state2.ξ20
    state1.ξ11[1:end] = state2.ξ11[1:end]
    for j in 1:timeindex
        state1.ξ02[timeindex,j] = state2.ξ02[timeindex,j]
    end
    for j in timeindex+1:length(state1.ξ11)
        state1.ξ02[j,timeindex] = state2.ξ02[j,timeindex]
    end
end

function set_zero!(state1::ZeroPhoton,timeindex::Int)
    state1.ξ0 = 0
end

function set_zero!(state1::OnePhoton,timeindex::Int)
    state1.ξ10 = 0
    state1.ξ01[timeindex] = 0
end

function set_zero!(state1::TwoPhoton,timeindex::Int)
    state1.ξ20 = 0
    state1.ξ11 .= 0
    for j in 1:timeindex
        state1.ξ02[timeindex,j] = 0
    end
    for j in timeindex+1:length(state1.ξ11)
        state1.ξ02[j,timeindex] = 0
    end
end



function mul!(state1::ZeroPhoton,a::ComplexF64,timeindex::Int)
    state1.ξ0 = state1.ξ0*a
end

function mul!(state1::OnePhoton,a::ComplexF64,timeindex::Int)
    state1.ξ10 = a*state1.ξ10
    state1.ξ01[timeindex] = state1.ξ01[timeindex]*a
end

function mul!(state1::TwoPhoton,a::ComplexF64,timeindex::Int)
    state1.ξ20 = state1.ξ20*a
    state1.ξ11 = state1.ξ11 .* a
    for j in 1:timeindex
        state1.ξ02[timeindex,j] = state1.ξ02[timeindex,j]*a
    end
    for j in timeindex+1:length(state1.ξ11)
        state1.ξ02[j,timeindex] = state1.ξ02[j,timeindex]*a
    end
end




function adw!(state::ZeroPhoton,diff::ZeroPhoton,timeindex::Int;g=1)
    "NOTHING"
end

function adw!(state::OnePhoton,diff::OnePhoton,timeindex;g=1)
    diff.ξ10 = diff.ξ10 + complex(g*state.ξ01[timeindex])
end

function adw!(state::TwoPhoton,diff::TwoPhoton,timeindex::Int;g=1)
    diff.ξ20 = diff.ξ20 + sqrt(2)*complex(g*state.ξ11[timeindex])
    for j in 1:timeindex-1
        diff.ξ11[j] =diff.ξ11[j] + complex(g*state.ξ02[timeindex,j])
    end
    for j in timeindex+1:length(state.ξ02[timeindex,:])
        diff.ξ11[j] =diff.ξ11[j]+complex(g*state.ξ02[j,timeindex])
    end
    diff.ξ11[timeindex] =diff.ξ11[timeindex]+complex(sqrt(2)*g*state.ξ02[timeindex,timeindex])
end

function wda!(state::ZeroPhoton,diff::ZeroPhoton,timeindex::Int;g=1)
    "NOTHING"
end

function wda!(state::OnePhoton,diff::OnePhoton,timeindex::Int;g=1)
    diff.ξ01[timeindex] = diff.ξ01[timeindex] + complex(g*state.ξ10)
end

function wda!(state::TwoPhoton,diff::TwoPhoton,timeindex::Int;g=1)
    diff.ξ11[timeindex] = diff.ξ11[timeindex] + sqrt(2)*complex(g*state.ξ20)
    for j in 1:timeindex-1
        diff.ξ02[timeindex,j] = diff.ξ02[timeindex,j]+complex(g*state.ξ11[j])
    end
    for j in timeindex+1:length(state.ξ11)
        diff.ξ02[j,timeindex] = diff.ξ02[j,timeindex] +complex(g*state.ξ11[j])
    end
    diff.ξ02[timeindex,timeindex] = diff.ξ02[timeindex,timeindex] + sqrt(2)*complex(g*state.ξ11[timeindex])
end


function annihilation!(state::TwoPhoton,diff::OnePhoton,timeindex::Int;g=1)
    diff.ξ10 = diff.ξ10 + complex(sqrt(2)*g*state.ξ20)
    diff.ξ01 = diff.ξ01 + complex(g*state.ξ11)  
end


function annihilation!(state::OnePhoton,diff::ZeroPhoton,timeindex::Int;g=1)
    diff.ξ0 = diff.ξ0 + complex(g*state.ξ10) 
end

function creation!(state::OnePhoton,diff::TwoPhoton,timeindex::Int;g=1)
    diff.ξ20 = diff.ξ20 + complex(sqrt(2)*g*state.ξ10)
    diff.ξ11 = diff.ξ11 + complex(g*state.ξ01)
end

function creation!(state::ZeroPhoton,diff::OnePhoton,timeindex::Int;g=1)
    diff.ξ10 = diff.ξ10 + complex(g*state.ξ0)
end

function number!(state::ZeroPhoton,diff::ZeroPhoton,timeindex::Int;g=1)
    "nothing"
end

function number!(state::OnePhoton,diff::OnePhoton,timeindex::Int;g=1)
    diff.ξ10 = diff.ξ10 + state.ξ10*complex(g) 
end
function number!(state::TwoPhoton,diff::TwoPhoton,timeindex::Int;g=1)
    diff.ξ20 = diff.ξ20 + 2*state.ξ20*complex(g)
    diff.ξ11 = diff.ξ11 .+ (state.ξ11 .* complex(g))
end

function annihilation!(state::CWState;g=1)
    for i in 1:length(state.data)-1
        annihilation!(state.data[i+1],state.diff[i],state.timeindex,g=g)
    end
end

function creation!(state::CWState;g=1)
    for i in 1:length(state.data)-1
        creation!(state.data[i],state.diff[i+1],state.timeindex,g=g)
    end
end

function wda!(state::CWState;g=1)
    for i in 1:length(state.data)
        wda!(state.data[i],state.diff[i],state.timeindex,g=g)
    end
end

function adw!(state::CWState;g=1)
    for i in 1:length(state.data)
        adw!(state.data[i],state.diff[i],state.timeindex,g=g)
    end
end

function number!(state::CWState;g=1)
    for i in 1:length(state.data)
        number!(state.data[i],state.diff[i],state.timeindex,g=g)
    end
end


function add_diff!(state1::CWState,state2::CWState)
    if state1.timeindex != state2.timeindex
        print("Carefull, timeindeces not same")
    end
    for i in 1:length(state1.data)
        add!(state1.data[i],state2.diff[i],state1.timeindex)
    end
end

function zero_diff!(state::CWState)
    for i in 1:length(state.data)
        zero_diff!(state.diff[i],state.timeindex)
    end
end

function add_diff!(state1::CWState)
    for i in 1:length(state1.data)
        add!(state1.data[i],state1.diff[i],state1.timeindex)
    end
end

function set_equal!(diff::CWState,state1::CWState)
    for i in 1:length(state1.data)
        set_equal!(diff.data[i],state1.data[i],diff.timeindex)
    end
    diff.timeindex = state1.timeindex 
end

function set_equal_diff!(state::CWState)
    for i in 1:length(state.data)
        set_equal!(state.data[i],state.diff[i],state.timeindex)
        zero_diff!(state.diff[i],state.timeindex)
    end
end

function timeindex_update!(state::CWState)
    state.timeindex = state.timeindex+1
end

function save_state!(state::CWState)
    if state.basis.N_cutoff == 0
        state.tmp_data[1] = state.data[1].ξ0
    elseif state.basis.N_cutoff == 1
        state.tmp_data[1] = state.data[1].ξ0
        state.tmp_data[2] = state.data[2].ξ10
        state.tmp_data[3:end] = state.data[2].ξ01 
    elseif state.basis.N_cutoff == 2
        N = length(state.basis.times)
        state.tmp_data[1] = state.data[1].ξ0
        state.tmp_data[2] = state.data[2].ξ10
        state.tmp_data[3:3+N-1] = state.data[2].ξ01
        state.tmp_data[3+N] = state.data[3].ξ20
        state.tmp_data[4+N:4+2N-1] = state.data[3].ξ11
        state.tmp_data[4+2N:4+2N+state.timeindex-1] = state.data[3].ξ02[state.timeindex,1:state.timeindex]
        state.tmp_data[4+2N+state.timeindex:end] = state.data[3].ξ02[state.timeindex+1:end,state.timeindex]
    end
end


function load_state!(state::CWState)
    if state.basis.N_cutoff == 0
        state.data[1].ξ0 = state.tmp_data[1]
    elseif state.basis.N_cutoff == 1
        state.data[1].ξ0 = state.tmp_data[1]
        state.data[2].ξ10 = state.tmp_data[2]
        state.data[2].ξ01 = state.tmp_data[3:end] 
    elseif state.basis.N_cutoff == 2
        N = length(state.basis.times)
        state.data[1].ξ0 = state.tmp_data[1]
        state.data[2].ξ10 = state.tmp_data[2]
        state.data[2].ξ01 = state.tmp_data[3:3+N-1]
        state.data[3].ξ20 = state.tmp_data[3+N]
        state.data[3].ξ11 = state.tmp_data[4+N:4+2N-1]
        state.data[3].ξ02[state.timeindex,1:state.timeindex] = state.tmp_data[4+2N:4+2N+state.timeindex-1]
        state.data[3].ξ02[state.timeindex+1:end,state.timeindex] = state.tmp_data[4+2N+state.timeindex:end]
    end
end


function save_diff!(state::CWState)
    if state.basis.N_cutoff == 0
        state.tmp_diff[1] = state.diff[1].ξ0
    elseif state.basis.N_cutoff == 1
        state.tmp_diff[1] = state.diff[1].ξ0
        state.tmp_diff[2] = state.diff[2].ξ10
        state.tmp_diff[3:end] = state.diff[2].ξ01 
    elseif state.basis.N_cutoff == 2
        N = length(state.basis.times)
        state.tmp_diff[1] = state.diff[1].ξ0
        state.tmp_diff[2] = state.diff[2].ξ10
        state.tmp_diff[3:3+N-1] = state.diff[2].ξ01
        state.tmp_diff[3+N] = state.diff[3].ξ20
        state.tmp_diff[4+N:4+2N-1] = state.diff[3].ξ11
        state.tmp_diff[4+2N:4+2N+state.timeindex-1] = state.diff[3].ξ02[state.timeindex,1:state.timeindex]
        state.tmp_diff[4+2N+state.timeindex:end] = state.diff[3].ξ02[state.timeindex+1:end,state.timeindex]
    end
end


function load_diff!(state::CWState)
    if state.basis.N_cutoff == 0
        state.diff[1].ξ0 = state.tmp_diff[1] + state.diff[1].ξ0
    elseif state.basis.N_cutoff == 1
        state.diff[1].ξ0 = state.tmp_diff[1] + state.diff[1].ξ0
        state.diff[2].ξ10 = state.tmp_diff[2] + state.diff[2].ξ10
        state.diff[2].ξ01 = state.tmp_diff[3:end] + state.diff[2].ξ01
    elseif state.basis.N_cutoff == 2
        N = length(state.basis.times)
        state.diff[1].ξ0 = state.tmp_diff[1] + state.diff[1].ξ0
        state.diff[2].ξ10 = state.tmp_diff[2] + state.diff[2].ξ10
        state.diff[2].ξ01 = state.tmp_diff[3:3+N-1] + state.diff[2].ξ01
        state.diff[3].ξ20 = state.tmp_diff[3+N] + state.diff[3].ξ20
        state.diff[3].ξ11 = state.tmp_diff[4+N:4+2N-1] + state.diff[3].ξ11
        state.diff[3].ξ02[state.timeindex,1:state.timeindex] = state.tmp_diff[4+2N:4+2N+state.timeindex-1] + state.diff[3].ξ02[state.timeindex,1:state.timeindex]
        state.diff[3].ξ02[state.timeindex+1:end,state.timeindex] = state.tmp_diff[4+2N+state.timeindex:end] + state.diff[3].ξ02[state.timeindex+1:end,state.timeindex]
    end
end


