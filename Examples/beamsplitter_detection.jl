using CavityWaveguide
using QuantumOptics

times = 0:0.1:20
bw = WaveguideBasis(1,times)
wa = destroy(bw)
wb = destroy(bw)

detector_plus = Detector(wa/sqrt(2),wb/sqrt(2))
detector_minus = Detector(wa/sqrt(2),-wb/sqrt(2))

ξfun(t,σ,t0) = complex(sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2))

waveguide_a = onephoton(bw,ξfun,1,10,times)
waveguide_b = onephoton(bw,ξfun,1,10,times)
ψ_total = LazyTensorKet(waveguide_a,waveguide_b)

p_plus_click = detect_single_click(ψ_total,detector_plus)
p_minus_click = detect_single_click(ψ_total,detector_minus)

println("Probability of having only one click in detector plus: $p_plus_click")
println("Probability of having only one click in detector plus: $p_minus_click")

p_plus_plus_click = detect_double_click(ψ_total,detector_plus,detector_plus)
p_minus_minus_click = detect_double_click(ψ_total,detector_minus,detector_minus)
p_plus_minus_click = detect_double_click(ψ_total,detector_plus,detector_minus)
p_minus_plus_click = detect_double_click(ψ_total,detector_minus,detector_plus)

println("Probability of having two clicks in detector plus: $p_plus_plus_click")
println("Probability of having two clicks in detector minus: $p_minus_minus_click")
println("Probability of having one click in detector plus and one in detector minus: $(p_plus_minus_click+p_minus_plus_click)")
