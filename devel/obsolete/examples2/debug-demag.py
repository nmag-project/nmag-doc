#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual

import ocaml #need this as long as we use ocaml.probe_field

from nmag2.si_units import SI

#produce more output for now
#nmag.set_global_log_level('debug')
#nmag.set_log_level('debug')

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
sim.timestepper_tuning_params=[1e-6,1e-6,2,300]

sim.defregion("Py", nm.ellipsoid([3.0,3.0,3.0]), mag_mat=mat_Py)

def initial_magnetization_0(coords,mag_type):
    return [0.8*math.sin(coords[0]/1.0)*1e6,
            0.8*math.cos(coords[0]/1.0)*1e6,
            0.6*1e6
            ]

def initial_magnetization(coords,mag_type):
    return [x*(-1.0) for x in coords] # radially outward-pointing. Will be normalized anyway.

# surface contribution is about +156
# bulk is about -109
# Difference: +46
#
# Analytical expectation: div e_r = 2/r, hence volume integral = 4pi R^2, gives 113.09

sim.generate_mesh(([-5.0,-5.0,-5.0],[5.0,5.0,5.0]), # bounding box
                  a0=1.0,
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

sim.advance_time(intensive_param_by_name,1e-6)

print "### XXX H_demag: ",ocaml.probe_field(sim.fields["H_demag"],"H_demag",[0.0,0.0,0.0])

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
