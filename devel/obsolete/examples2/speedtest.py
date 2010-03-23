#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere
# Speed test

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual

import ocaml #need this as long as we use ocaml.probe_field

from nmag2.si_units import SI

time_init=time.time()

#produce more output for now
#nmag.set_global_log_level('debug')
#nmag.set_log_level('debug')

intensive_param_by_name={"H_x":0.1,"H_y":0.0,"H_z":0.0}
# very slightly pulling in x-direction

mat_Py = nmag.MagMaterial("Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12,"J/m"),
                          )

sim=nmag.SimulationContext("sphere")

# sim.timestepper_tuning_params=[1e-6,1e-6,2,300] # These are the defaults...
sim.timestepper_tuning_params=[1e-6,1e-6,2,300]

sim.defregion("Py", nm.ellipsoid([3.0,3.0,3.0]), mag_mat=mat_Py)

def initial_magnetization(coords,mag_type):
    return [0.8*math.sin(coords[0]/1.0)*1e6,
            0.8*math.cos(coords[0]/1.0)*1e6,
            0.6*1e6
            ]

sim.generate_mesh(([-5.0,-5.0,-5.0],[5.0,5.0,5.0]), # bounding box
                  a0=0.7,
                  max_steps=450,
                  cache_name="sphere"
                  )

#sim.set_magnetization([1.0,0.0,0.0])
#sim.set_magnetization([0.0,1.0,0.0])
#sim.set_magnetization([0.0,0.0,1.0])
#sim.set_magnetization([0.0,1.0,1.0])
# ^All these yield a magnetic charge imbalance of 0.0


sim.set_magnetization(initial_magnetization)

time_start=time.time()

sim.advance_time(intensive_param_by_name,2.0)
#sim.advance_time(intensive_param_by_name,10.0)

print "H_demag at origin: ",ocaml.probe_field(sim.fields["H_demag"],"H_demag",[0.0,0.0,0.0])


time_end=time.time()

print "SIM INIT TIME: ",time_start-time_init," sec."
print "SIM RUN TIME: ",time_end-time_start," sec."
print "PETSC: ",ocaml.petsc_cpu_cycles("report")
print "PETSC: ",ocaml.petsc_cpu_cycles("reset")
sys.exit()

print
"""

nmag, basic run, mesh pre-generated, cached:

DDD CPU cvode_advance total: 213505850182.000000 CPU cycles
SIM INIT TIME:  106.649693012  sec.
SIM RUN TIME:  226.640647888  sec.

Getting the computation of rho and phi right:

SIM INIT TIME:  106.865993023  sec.
SIM RUN TIME:  164.981091022  sec.

Ad influence of stack/GC parameters:

CAMLRUNPARAM="s=1000k,i=4000k" ../bin/nsim speedtest.py

DDD CPU cvode_advance total: 150980343400.000000 CPU cycles
SIM INIT TIME:  105.664279938  sec.
SIM RUN TIME:  158.250349045  sec.


Going back to more primitive Jacobi computation:

DDD CPU cvode_advance total: 189762587372.000000 CPU cycles
SIM INIT TIME:  107.145243883  sec.
SIM RUN TIME:  200.573558092  sec.

DDD CPU cvode_advance total: 198934409488.000000 CPU cycles
SIM INIT TIME:  101.541635036  sec.
SIM RUN TIME:  209.417989016  sec.

~~~~~ Everything properly optimized again:

DDD CPU cvode_advance total: 137655170811.000000 CPU cycles
SIM INIT TIME:  108.317420959  sec.
SIM RUN TIME:  148.439527035  sec.

DDD CPU cvode_advance total: 138210050850.000000 CPU cycles
SIM INIT TIME:  107.049017906  sec.
SIM RUN TIME:  148.985625982  sec.

Changing FEM symbolic algebra to use caching throughout:

DDD CPU cvode_advance total: 145667462481.000000 CPU cycles
SIM INIT TIME:  101.825253963  sec.
SIM RUN TIME:  156.203269005  sec.

...hence, some little advantage from that (and therefore, I leave it as it is),
but the vast majority of our setup time is being burnt somewhere else!

Note ad Jacobian computation times:

MatCreate:
DDD CPU fun_jacobi total: 56321946966.000000 CPU cycles

MatCreateSeqAIJ: (pre_allocation=50):
DDD CPU fun_jacobi total:  2995325250.000000 CPU cycles

Speedup factor: 18.8

========

Setup time:

DDD diffop_prematrix T=13.828069
Operator:
#<Differential Operator 
 -   20.69014 * <d/dx2 H_exch_Py[0] || d/dx2 m_Py[0]>
 -   20.69014 * <d/dx1 H_exch_Py[0] || d/dx1 m_Py[0]>
 -   20.69014 * <d/dx0 H_exch_Py[0] || d/dx0 m_Py[0]>
 -   20.69014 * <d/dx2 H_exch_Py[1] || d/dx2 m_Py[1]>
 -   20.69014 * <d/dx1 H_exch_Py[1] || d/dx1 m_Py[1]>
 -   20.69014 * <d/dx0 H_exch_Py[1] || d/dx0 m_Py[1]>
 -   20.69014 * <d/dx2 H_exch_Py[2] || d/dx2 m_Py[2]>
 -   20.69014 * <d/dx1 H_exch_Py[2] || d/dx1 m_Py[2]>
 -   20.69014 * <d/dx0 H_exch_Py[2] || d/dx0 m_Py[2]>

>


DDD diffop_prematrix T=4.275483
Operator:
#<Differential Operator 
 +    1.00000 * <H_exch_Py[0] || H_exch_Py[0]>
 +    1.00000 * <H_exch_Py[1] || H_exch_Py[1]>
 +    1.00000 * <H_exch_Py[2] || H_exch_Py[2]>

>

DDD diffop_prematrix T=5.658869
Operator:
#<Differential Operator 
 +    1.00000 * <rho || d/dx0 m_Py[0]>
 +    1.00000 * <rho || d/dx1 m_Py[1]>
 +    1.00000 * <rho || d/dx2 m_Py[2]>

>

DDD diffop_prematrix T=4.448848
Operator:
#<Differential Operator 
 -    1.00000 * <d/dx2 rho || d/dx2 phi>
 -    1.00000 * <d/dx1 rho || d/dx1 phi>
 -    1.00000 * <d/dx0 rho || d/dx0 phi>

>

DDD diffop_prematrix T=6.078495
Operator:
#<Differential Operator 
 -    1.00000 * <H_demag[0] || d/dx0 phi>
 -    1.00000 * <H_demag[1] || d/dx1 phi>
 -    1.00000 * <H_demag[2] || d/dx2 phi>

>

DDD diffop_prematrix T=3.748707
Operator:
#<Differential Operator 
 +    1.00000 * <H_demag[0] || H_demag[0]>
 +    1.00000 * <H_demag[1] || H_demag[1]>
 +    1.00000 * <H_demag[2] || H_demag[2]>

>

Total diffop_prematrix operator setup time: 37.768 seconds
(hence, approx. 40%)


High precision results:

H_demag at origin:  [('H_demag', [-0.025628112201029212, -0.14303234266097389, -0.29270063950231923])]

Reducing KSP precision to 1%:

DDD CPU cvode_advance total: 107768410388.000000 CPU cycles
H_demag at origin:  [('H_demag', [-0.030356189679353454, -0.13806744317357139, -0.29121662009153099])]
SIM INIT TIME:  58.1213579178  sec.
SIM RUN TIME:  117.946659088  sec.
PETSC:  (('Matrix Solver', 66083744840.0), ('Matrix*Vector', 13831156955.0), ('Vector', 168225857.0))

... speed increase by another ~22% - presumably not worth it?

"""

