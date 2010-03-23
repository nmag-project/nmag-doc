

#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere

import nmag2 as nmag, nmesh as nmesh, math, time
from nmag2.si_units import SI
import nfem
import nfem.visual

import ocaml #need this as long as we use ocaml.probe_field



intensive_param_by_name={"H_x":0.2,"H_y":0.0,"H_z":0.0}
# very slightly pulling in x-direction

def test_anisotropy_energy(v): # direction is unit-normalized!
    print "DDD anisotropy dir=",v
    result=(1.0*v[0]*v[0]+2.0*v[1]*v[1]+3.5*v[2]*v[2])
    return result

mat_Py = nmag.MagMaterial("Py",
                          Ms=1.0,
                          J=13.0,
                          anisotropy_order=2,
                          anisotropy=test_anisotropy_energy
                          )

mat_Fe = nmag.MagMaterial("Fe",
                          Ms=SI(1.7e6,"A/m"),
                          J=SI(2.07e-11,"J/m")
                          )

sim=nmag.SimulationContext("sphere")

# sim.timestepper_tuning_params=[1e-6,1e-6,2,300] # These are the defaults...
sim.timestepper_tuning_params=[1e-6,1e-6,4,300]

sim.set_magnetization([0.0,1.0,0.0])
# ^ just to show we can place set_magnetization() where we want...

sim.defregion("Py", nmesh.ellipsoid([2.0,2.0,2.0],transform=[('shift',[5,0,0])]), mag_mat=mat_Py)
sim.defregion("Fe", nmesh.ellipsoid([2.0,2.0,2.0]), mag_mat=mat_Fe)


def initial_magnetization(coords,mag_type):
    if mag_type == "m_Fe":
        return [1,0,0]
    else:
        return [0.8*math.sin(coords[0]/3.0)*1e6,
                0.8*math.cos(coords[0]/3.0)*1e6,
                0.6*1e6
                ]

sim.generate_mesh(([-2.5,-2.5,-2.5],[7.5,2.5,2.5]), # bounding box
                  a0=1.0,
                  max_steps=1200,
                  cache_name="sphere3"
                  )

sim.mesh.save('twomaterials.nmesh')

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

nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-0.vtk',mesh=sim.mesh)

time_start=time.time()


for i in range(1,2):
    target_time=sim.advance_time(intensive_param_by_name,0.05*i)
    #write SI data into default file (results.ndt)
    sim.save_data_table()
    #write SU data into special file 
    sim.save_data_table_su(filename='results_su.ndt')
    nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-%d.vtk' % i,mesh=sim.mesh)
    print "Field at origin (T=%f):"%target_time,sim.get_m([0.0,0.0,0.0])
time_end=time.time()

print "SIM RUN TIME: ",time_end-time_start," sec."

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


