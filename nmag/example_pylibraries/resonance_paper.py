# Shortened/prettifyed version of the script for the nmag-par paper

import nmag, math, scipy.optimize
from nmag import SI, every

Py = nmag.MagMaterial(name = 'Py',
                      Ms = SI(0.86e6, 'A/m'),
                      llg_damping=0.02,
                      exchange_coupling=
                      SI(13.0e-12, 'J/m'))

A_mag = SI(0.001e6, "A/m") # amplitude
sim = nmag.Simulation()
sim.load_mesh('rod.nmesh.h5',[('rod', Py)],
              unit_length=SI(1e-9,'m'))

def simulate_resonance(freq,
                       amplitude=A_mag,
                       dt=SI(5*1e-12,"s"),
                       time_steps=800
                       ):
  sim.set_m([0,0,1])
  angle_step=(2*math.pi*freq*dt).value
  result=[]
  for i in range(time_steps):
    h_ext=[amplitude.value*
           math.sin(i*angle_step),0,0]
    sim.set_H_ext(h_ext,SI(1,"A/m"))
    sim.advance_time(i*dt)
    a=sim.get_subfield_average_siv("M","Py")
    result.append((i*dt,h_ext,a))
  return result
                       
def amplitude(f_GHz):
  res=simulate_resonance(SI(f_GHz*1e9,"1/s"))
  osc=[[r[0].value,r[2][0],r[2][1],r[2][2]]
       for r in res[400:]] # skip lead_in
  (f,params)=nmag.fit_oscillation(osc)
  a=math.sqrt(sum([p[1]*p[1] for p in params]))
  return -a # we minimize -(amplitude)!

freq=scipy.optimize.fmin(amplitude,
                         [1.5], # start: 1.5 GHz
                         xtol=0.05,
                         ftol=0.05)
print "Resonance frequency: %.1f GHz"%freq
