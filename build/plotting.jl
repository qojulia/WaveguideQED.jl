"""
    plot_twophoton!(ax,twophotonstate::TwophotonView,times)
    plot_twophoton!(ax,state::Ket,times)

Plots the twophoton state in the given ax. If state is a `Ket` [`TwoPhotonView`](@ref) is called to extract twophotonstate. 
Returns ax.contour object.
"""
function plot_twophoton!(ax,twophotonstate::TwoPhotonView,times)
    xgrid = repeat(times',length(times),1)
    ygrid = repeat(times,1,length(times))
    cnt1= ax.contourf(xgrid,ygrid,twophotonstate.*conj(twophotonstate),100)
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    cnt1
end
function plot_twophoton!(ax,state::Ket,times)
    twophotonstate = TwoPhotonView(state)
    xgrid = repeat(times',length(times),1)
    ygrid = repeat(times,1,length(times))
    cnt1= ax.contourf(xgrid,ygrid,twophotonstate.*conj(twophotonstate),100)
    for c in cnt1.collections
        c.set_edgecolor("face")
    end
    ax.set_aspect("equal", "box")
    cnt1
end