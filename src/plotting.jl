"""
    plot_twophoton!(ax,twophotonstate::TwophotonView,times;kwargs...)
    plot_twophoton!(ax,twophotonstate::TwoWaveguideView,times;kwargs...)
    plot_twophoton!(ax,state::Ket,times;kwargs...)
    plot_twophoton!(ax,twophotonstate,times;kwargs...)

Plots the twophoton state in the given ax. 
# Arguments
- ax of type PyObject <AxesSubplot: > from `PyPlot`
- State to be plotted twophotonstate or state. If state is a `Ket` [`TwoPhotonView`](@ref) is called to extract twophotonstate. Otherwise twophotonstate should be AbstractArray of dimensions (length(times),length(times)).

#Return
ax.contour object with the plotted state.
"""
function plot_twophoton!(ax,twophotonstate::TwoPhotonView,times;kwargs...)
    xgrid = repeat(times',length(times),1)
    ygrid = repeat(times,1,length(times))
    cnt1= ax.contourf(xgrid,ygrid,abs.(twophotonstate).^2,100;kwargs...)
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    cnt1
end
function plot_twophoton!(ax,twophotonstate::TwoWaveguideView,times;kwargs...)
    xgrid = repeat(times',length(times),1)
    ygrid = repeat(times,1,length(times))
    cnt1= ax.contourf(xgrid,ygrid,abs.(twophotonstate).^2,100;kwargs...)
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    cnt1
end
function plot_twophoton!(ax,state::Ket,times;kwargs...)
    twophotonstate = TwoPhotonView(state)
    xgrid = repeat(times',length(times),1)
    ygrid = repeat(times,1,length(times))
    cnt1= ax.contourf(xgrid,ygrid,abs.(twophotonstate).^2,100;kwargs...)
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    cnt1
end
function plot_twophoton!(ax,twophotonstate,times;kwargs...)
    xgrid = repeat(times',length(times),1)
    ygrid = repeat(times,1,length(times))
    cnt1= ax.contourf(xgrid,ygrid,abs.(twophotonstate).^2,100;kwargs...)
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    cnt1
end