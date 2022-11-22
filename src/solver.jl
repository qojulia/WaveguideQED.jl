export cwsolver!

#The timestep algorithm using a RK4 method. 
function timestep!(state::CWState,k1::CWState,k2::CWState,k3::CWState,k4::CWState,U::CWOperator,dt,timeindex)
    #Compute k1 derivative
    mul!(k1,U)
    mul_diff!(k1,(dt/2))
    
    #Add the k1 derivative to k2 and compute derivative
    add_diff!(k2,k1)
    mul!(k2,U)
    mul_diff!(k2,(dt/2))
    
    #Add the k2 derivative to k3 and compute derivative
    add_diff!(k3,k2)
    mul!(k3,U)
    mul_diff!(k3,dt)
    
    #Add the k3 derivative to k4 and compute derivative
    add_diff!(k4,k3)
    mul!(k4,U)
    
    #Multiple all derivatives with factors for adding to final state
    mul_diff!(k1,(1/3))
    mul_diff!(k2,(2/3))
    mul_diff!(k3,(1/3))
    mul_diff!(k4,(dt/6))

    #Add all derivatives to the output state
    add_diff!(state,k1)
    zero_diff!(k1)
    add_diff!(state,k2)
    zero_diff!(k2)
    add_diff!(state,k3)
    zero_diff!(k3)
    add_diff!(state,k4)
    zero_diff!(k4)

    #Update k1,k2,k3, and k4 to match output state for next time step
    timeindex_update!(state)
    #state.timeindex = timeindex+1
    set_equal!(k1,state)
    set_equal!(k2,state)
    set_equal!(k3,state)
    set_equal!(k4,state)
end


function cwsolver!(H::T,ξin::CWState) where T<:CWOperator
    dt = ξin.basis.times[2] - ξin.basis.times[1] 
    k1 =  copy(ξin)
    k2 =  copy(ξin)
    k3 =  copy(ξin)
    k4 =  copy(ξin)
    # Cavity time evolution operator (first order)
    U = -im * H
    for i in 1:length(ξin.basis.times)
        timestep!(ξin,k1,k2,k3,k4,U,dt,i)
    end
end


