using Test
using WaveguideQED
using QuantumOptics
include("helper_functions.jl")

param = BarretKokParameters()
param.times = 0:0.1:30
param.γ = 0
sys = prep_fast(param)
bc = sys.bc
bw = sys.bw
be = sys.be

ξfun(t,σ,t0) = complex(sqrt(2/σ)* (log(2)/pi)^(1/4)*exp(-2*log(2)*(t-t0)^2/σ^2))
psi_a_1 = onephoton(bw,ξfun,param.times,1,15) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,1)
psi_b_1 = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,2)
psi_a_2 = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,2)
psi_b_2 = onephoton(bw,ξfun,param.times,1,15) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,1)
proj_up = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,2)
proj_down = zerophoton(bw) ⊗ fockstate(bc,0) ⊗ nlevelstate(be,1)

p1 = LazyTensorKet(proj_up,proj_down)
p2 = LazyTensorKet(proj_down,proj_up)

psi1 = LazyTensorKet(psi_a_1,psi_b_1)
psi2 = LazyTensorKet(psi_a_2,psi_b_2)

Detector_plus = Detector(sys.wa/sqrt(2),sys.wb/sqrt(2))
Detector_minus = Detector(sys.wa/sqrt(2),-sys.wb/sqrt(2))

_,p = detect_single_click(psi1,Detector_plus,p2)
@test isapprox(p,0.5)
_,p = detect_single_click(psi1,Detector_plus,p2)
@test isapprox(p,0.5)

_,p = detect_single_click(psi1,Detector_minus,p2)
@test isapprox(p,0.5)

_,p = detect_single_click(psi2,Detector_plus,p1)
@test isapprox(p,0.5)

_,p = detect_single_click(psi2,Detector_minus,p1)
@test isapprox(p,0.5)

_,p = detect_single_click(psi1,Detector_plus,p1)
@test isapprox(p+1,1)

_,p = detect_single_click(psi1,Detector_minus,p1)
@test isapprox(p+1,1)

_,p = detect_single_click(psi2,Detector_plus,p2)
@test isapprox(p+1,1)

_,p = detect_single_click(psi2,Detector_minus,p2)
@test isapprox(p+1,1)


psi_plus = (psi1+psi2)/sqrt(2)
psi_minus = (psi1-psi2)/sqrt(2)
projector_plus = (p1+p2)/sqrt(2)
projector_minus = (p1-p2)/sqrt(2)

_,p = detect_single_click(psi_plus,Detector_plus,projector_plus)
@test isapprox(p,0.5)

_,p = detect_single_click(psi_plus,Detector_minus,projector_minus)
@test isapprox(p,0.5)

_,p = detect_single_click(psi_plus,Detector_plus,projector_minus)
@test isapprox(p+1,1)

_,p = detect_single_click(psi_plus,Detector_minus,projector_plus)
@test isapprox(p+1,1)

_,p = detect_single_click(psi_minus,Detector_plus,projector_plus)
@test isapprox(p+1,1)

_,p = detect_single_click(psi_minus,Detector_minus,projector_minus)
@test isapprox(p+1,1)

_,p = detect_single_click(psi_minus,Detector_plus,projector_minus)
@test isapprox(p,0.5)

_,p = detect_single_click(psi_minus,Detector_minus,projector_plus)
@test isapprox(p,0.5)
