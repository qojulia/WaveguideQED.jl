export cwsolver!

#The timestep algorithm using a RK4 method. 
function timestep!(state::cwstate,k1::cwstate,k2::cwstate,k3::cwstate,k4::cwstate,H::CWOperator,dt,timeindex)
    # Cavity time evolution operator (first order)
    U = -im * H
    #Compute k1 derivatives
    U*k1
    mul_diff!(k1,(dt/2))
    
    add_diff!(k2,k1)
    U*k2
    mul_diff!(k2,(dt/2))
    
    add_diff!(k3,k2)
    U*k3
    mul_diff!(k3,dt)
    
    add_diff!(k4,k3)
    U*k4
    
    mul_diff!(k1,(1/3))
    mul_diff!(k2,(2/3))
    mul_diff!(k3,(1/3))
    mul_diff!(k4,(dt/6))

    add_diff!(state,k1)
    zero_diff!(k1)
    add_diff!(state,k2)
    zero_diff!(k2)
    add_diff!(state,k3)
    zero_diff!(k3)
    add_diff!(state,k4)
    zero_diff!(k4)
    state.timeindex = timeindex+1
    set_equal!(k1,state)
    set_equal!(k2,state)
    set_equal!(k3,state)
    set_equal!(k4,state)
end


function cwsolver!(H::T,ξin::cwstate) where T<:CWOperator
    dt = ξin.basis.times[2] - ξin.basis.times[1] 
    k1 =  deepcopy(ξin)
    k2 =  deepcopy(ξin)
    k3 =  deepcopy(ξin)
    k4 =  deepcopy(ξin)
    for i in 1:length(ξin.basis.times)
        timestep!(ξin,k1,k2,k3,k4,H,dt,i)
    end
end


