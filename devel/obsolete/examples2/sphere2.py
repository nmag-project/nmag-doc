#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual

intensive_param_by_name={"H_x":0.1,"H_y":0.0,"H_z":0.0}
# very slightly pulling in x-direction

mat_Py = nmag.MagMaterial("Py",Ms=1.0,J=13.0,  # A = exchange constant
                          alpha=0.5,
                          beta=0.2,
                          gamma=2.0,
                          )




sim=nmag.SimulationContext("sphere",
                           )
                           
sim.defregion("Py", nm.ellipsoid([3.0,3.0,3.0]), mag_mat=mat_Py)

def initial_magnetization(dof_name,coords):
    direction=dof_name[1][0]
    if direction==2:
        return 0.6
    elif direction==1:
        return 0.8*math.cos(coords[0]/3.0)
    else:
        return 0.8*math.sin(coords[0]/3.0)

sim.generate_mesh(([-5.0,-5.0,-5.0],[5.0,5.0,5.0]), # bounding box
                  initial_magnetization,
                  # extra meshing params
                  # a0=0.5,
                  a0=1.2,
                  max_steps=1000,
                  cache_name="sphere2"
                  )

# sim.set_magnetization(initial_magnetization)

print "Field at origin:",ocaml.probe_field(sim.fields["m"],"",[0.0,0.0,0.0])

# ...this shows that field_m was indeed set properly...

nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-0.vtk',mesh=sim.mesh)

for i in range(1,500):
    (target_time,target_field)=sim.advance_time(intensive_param_by_name,0.05*i)
    nfem.visual.fields2vtkfile(sim.fields.values(),'sphere-%d.vtk' % i,mesh=sim.mesh)

