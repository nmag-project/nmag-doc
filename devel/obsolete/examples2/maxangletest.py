#
# Copied from sphere3, and zeroed anisotropy and made initial magnetisation uniform.
#

#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual

import ocaml #need this as long as we use ocaml.probe_field

from nmag2.si_units import SI

import nmag2.si_units as si_units

#produce more output for now
nmag.set_global_log_level('debug')
nmag.set_log_level('debug')

intensive_param_by_name={"H_x":0.0,"H_y":0.0,"H_z":0.0}
# very slightly pulling in x-direction

def test_anisotropy_energy(v): # direction is unit-normalized!
    print "DDD anisotropy dir=",v
    result=(1.0*v[0]*v[0]+2.0*v[1]*v[1]+3.5*v[2]*v[2])
    return result*0.0

mat_Py = nmag.MagMaterial("Py",
                          Ms=SI(1e6,"A/m"),
                          J=SI(13.0e-12,"J/m"),
                          #Ms=1.,      
                          #J=13.,      
                          anisotropy_order=2,
                          anisotropy=test_anisotropy_energy
                          )

sim=nmag.SimulationContext("sphere")

# sim.timestepper_tuning_params=[1e-6,1e-6,2,300] # These are the defaults...
sim.timestepper_tuning_params=[1e-6,1e-6,4,300]

#sim.set_magnetization([0.0,1.0,0.0])
# ^ just to show we can place set_magnetization() where we want...



sim.defregion("Py", nm.ellipsoid([20.0,20.0,5.0]), mag_mat=mat_Py)

def initial_magnetization(coords,mag_type):
    #return [1,0,0]
    return [0.8*math.sin(coords[2]/5.0*6.28)*1e6,
            0.8*math.cos(coords[2]/5.0*6.28)*1e6,
            1e6
            ]

sim.set_magnetization(initial_magnetization)

sim.generate_mesh(([-25.0,-25.0,-25.0],[25.0,25.0,25.0]), # bounding box
                  a0=2.0,
                  max_steps=200,
                  cache_name="sphere3"
                  )

#nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-0.vtk',mesh=sim.mesh)

#sim.save_fields_vtk(filename='sphere-%d.vtk' % 0)

time_start=time.time()

for i in range(1,50):
    target_time=sim.advance_time(intensive_param_by_name,0.05*i)

    #write SI data into default file (results.ndt)
    sim.save_data_table()

    #write SU data into special file 
    sim.save_data_table_su(filename='results_su.ndt')

    sim.save_fields_vtk(filename='sphere-%d.vtk' % i)

    print "Field at origin (T=%f):"%target_time,sim.get_m([0.0,0.0,0.0])

time_end=time.time()

print "SIM RUN TIME: ",time_end-time_start," sec."

if False:
    f=sim._get_nfem_fields()
    print "Fields are",f

    print "probe magnetisation along line:"

    nm = 1e-9
    data = []
    for i in range(-10,11):
        x = float(i)/10*3*nm
        data.append(sim.get_M([x,0,0]))
    print data

    print "Demag field at [0,0,0] =", sim.get_H_demag([0,0,0])
    print "Exchange field at [0,0,0] =", sim.get_H_exchange([0,0,0])
    print "anisotropy field at [0,0,0] =", sim.get_H_anisotropy([0,0,0])

