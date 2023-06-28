# Waveguides with Time-Delay

```@example timedelay
times = 0:0.1:12 #Times for setting up waveguide basis
dt = times[2]-times[1]
be = FockBasis(1)
bw = WaveguideBasis(1,1,times)
sdw = create(be) |$\otimes$| destroy(bw)
wds = destroy(be) |$\otimes$| create(bw)
```

```
gamma = 1 #Decay rate of emitter
delay_time = 1 #In units of gamma
phi = pi #Mirror phase
``` 

Create delayed waveguide operator with delay keyword
```
sdw_delay = create(be) |$\otimes$| destroy(bw;delay=delay_time/dt)
wds_delay = destroy(be) |$\otimes$| create(bw;delay=delay_time/dt)
```

```
H = exp(i*phi)*sqrt(gamma/2/dt)*(sdw+wds)+sqrt(gamma/2/dt)*(sdw_delay+wds_delay)
```

Operator for emitter population

```
sd = create(be) |$\otimes$| identityoperator(bw)
s = destroy(be) |$\otimes$| identityoperator(bw)
n = ad*a
function ne_exp(time,psi)
    expect(n,psi)
end
```

```
psi_initial = fockstate(be,1) |$\otimes$| zerophoton(bw) #Initial state
times_sim = 0:0.1:10 #Simulation time has to be smaller than times due to delay. 
tout, ne_pi = waveguide_evolution(times_sim, psi_initial, H,fout=ne_exp)
```
