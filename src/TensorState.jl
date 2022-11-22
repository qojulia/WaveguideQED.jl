mutable struct CWTensorState
    ket::Ket
    ketbuffer::Ket
    cwstates::Vector{CWState}
end

function lazytensor(ket::Ket,cwstate::CWState)
    cwstates = Array{CWState}(undef,length(ket.data))
    for i in 1:length(ket.data)
        cwstates[i] = copy(cwstate)
    end
    return CWTensorState(ket,copy(ket),cwstates)
end

lazytensor(cwstate::CWState,ket::StateVector) = lazytensor(ket,cwstate)

copy(CWTensorState) = CWTensorState(copy(CWTensorState.ket),copy(CWTensorState.ketbuffer),copy(CWTensorState.cwstates))

function mul_diff!(state::CWTensorState,a::ComplexF64)
    state.ketbuffer = state.ketbuffer*a
    for i in 1:length(state.cwstates)
        mul_diff!(state.cwstate[i],a)
    end 
end

function add_diff!(state1::CWTensorState,state2::CWTensorState)
    state1.ket = state1.ket + state2.ketbuffer
    for i in 1:length(state1.cwstates)
        add_diff!(state1.cwstate[i],state2.cwstate[i])
    end 
end

function zero_diff!(state::CWTensorState)
    state.ketbuffer = 0*state.ketbuffer
    for i in 1:length(state.cwstates)
        zero_diff!(state.cwstate[i])
    end 
end

function set_equal!(diff::CWTensorState,state::CWTensorState)
    diff.ket = state.ket
    for i in 1:length(state.cwstates)
        set_equal!(diff.cwstate[i],state.cwstate[i])
    end
end

function timeindex_update!(state::CWTensorState)
    for i in 1:length(state.cwstates)
        timeindex_update!(state.cwstates[i])
    end
end

function get_output(state::CWTensorState)
    print("Output state occupation: $(state.ket.data[1])")    
    return  state.cwstates[1]
end