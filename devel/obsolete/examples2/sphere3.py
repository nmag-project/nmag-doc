#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual

import ocaml #need this as long as we use ocaml.probe_field

from nmag2.si_units import SI

#produce more output for now
nmag.set_global_log_level('debug')
nmag.set_log_level('debug')

intensive_param_by_name={"H_x":0.1,"H_y":0.0,"H_z":0.0}
# very slightly pulling in x-direction

def test_anisotropy_energy(v): # direction is unit-normalized!
    print "DDD anisotropy dir=",v
    result=(1.0*v[0]*v[0]+2.0*v[1]*v[1]+3.5*v[2]*v[2])
    return result*0.0

mat_Py = nmag.MagMaterial("Py",
                          Ms=SI(1e6,"A/m"),
                          J=SI(13.0e-12,"J/m"),
                          #Ms=1.,                #Matteo, very strange: if I use the SI units above,
                          #J=13.,                #then I get a wierd magsim-brain failure!
                          anisotropy_order=2,
                          anisotropy=test_anisotropy_energy
                          )

sim=nmag.SimulationContext("sphere")

# sim.timestepper_tuning_params=[1e-6,1e-6,2,300] # These are the defaults...
sim.timestepper_tuning_params=[1e-6,1e-6,4,300]

sim.set_magnetization([0.0,1.0,0.0])
# ^ just to show we can place set_magnetization() where we want...

sim.defregion("Py", nm.ellipsoid([4.0,4.0,5.0]), mag_mat=mat_Py)

def initial_magnetization(coords,mag_type):
    return [0.8*math.sin(coords[0]/1.0)*1e6,
            0.8*math.cos(coords[0]/1.0)*1e6,
            0.6*1e6
            ]

sim.generate_mesh(([-30.0,-30.0,-50.0],[30.0,30.0,50.0]), # bounding box
                  a0=2.0,
                  max_steps=400,
                  cache_name="sphere3"
                  )

#
#sim.set_magnetization(initial_magnetization)
#
#import ocaml
#print "Field at origin (initial):",ocaml.probe_field(sim.fields["m"],"",[0.0,0.0,0.0])
#
#for i in range(1,3):
#    target_time=sim.advance_time(intensive_param_by_name,0.05*i)
#    print "Field at origin (T=%f):"%target_time,ocaml.probe_field(sim.fields["m"],"",[0.0,0.0,0.0])

# Re-set magnetization:

sim.set_magnetization(initial_magnetization)

time_start=time.time()

for i in range(1,100):
    target_time=sim.advance_time(intensive_param_by_name,0.1*i)

    #write SI data into default file (results.ndt)
    sim.save_data_table()

    #write SU data into special file 
    sim.save_data_table_su(filename='results_su.ndt')

    nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-%d.vtk' % i,mesh=sim.mesh)
    print "Field at origin (T=%f):"%target_time,sim.get_m([0.0,0.0,0.0])
time_end=time.time()

print "SIM RUN TIME: ",time_end-time_start," sec."
