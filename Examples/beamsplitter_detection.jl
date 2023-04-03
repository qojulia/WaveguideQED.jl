using WaveguideQED
using QuantumOptics

times = 0:0.1:20
bw = WaveguideBasis(1,times)
wa = destroy(bw)
wb = destroy(bw)

detector_plus = Detector(wa/sqrt(2),wb/sqrt(2))
detector_minus = Detector(wa/sqrt(2),-wb/sqrt(2))

ξfun(t,σ,t0) = complex(sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2))

waveguide_a = onephoton(bw,ξfun,times,1,10)
waveguide_b = onephoton(bw,ξfun,times,1,10)
ψ_total = LazyTensorKet(waveguide_a,waveguide_b)

p_plus_click = detector_plus*ψ_total
p_minus_click = detector_minus*ψ_total

println("Probability of having only one click in detector plus: $p_plus_click")
println("Probability of having only one click in detector plus: $p_minus_click")

p_plus_plus_click = detector_plus*detector_plus*ψ_total
p_minus_minus_click = detector_minus*detector_minus*ψ_total
p_plus_minus_click = detector_minus*detector_plus*ψ_total
p_minus_plus_click = detector_plus*detector_minus*ψ_total

println("Probability of having two clicks in detector plus: $p_plus_plus_click")
println("Probability of having two clicks in detector minus: $p_minus_minus_click")
println("Probability of having one click in detector plus and one in detector minus: $(p_plus_minus_click+p_minus_plus_click)")

waveguide_a = onephoton(bw,ξfun,times,1,5)
waveguide_b = onephoton(bw,ξfun,times,1,15)
ψ_total = LazyTensorKet(waveguide_a,waveguide_b)
p_plus_plus_click = detector_plus*detector_plus*ψ_total
p_minus_minus_click = detector_minus*detector_minus*ψ_total
p_plus_minus_click = detector_minus*detector_plus*ψ_total
p_minus_plus_click = detector_plus*detector_minus*ψ_total
