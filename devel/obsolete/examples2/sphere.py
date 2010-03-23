#
# (C) 2006 Dr. Thomas Fischbacher
# Relaxation of the homogeneously magnetized sphere

import nmag1 as mag, nmesh as nm, sys, math, time
import nfem # for debugging

#mag.set_intensive_parameters(["H_x","H_y","H_z"], # may also set T,p,E, etc.
#                             external_magnetization=["H_x","H_y","H_z"])
# ^ XXX not used yet, as we are going to use the experimental cvode integrator...

mag.set_intensive_parameters([])


# mag.set_default_order(1) # default anyway...

mag.set_features({"demag":True,"exchange":True,"timestep":True}) # XXX change usage!

mat_Py = mag.MagMaterial("Py",Ms=1.0,A=13.0,  # A = exchange constant
                         )

mag.defregion("Py", nm.ellipsoid([3.0,3.0,3.0]), mag_mat=mat_Py)

mag.set_meshing_parameters(cache_name="exchange-spring-mesh",
                           bounding_box=([-4.0,-4.0,-4.0],[4.0,4.0,4.0]),
                           a0=1.0,
                           max_steps=600,
                           )

mag.create_mesh()

mag.set_magnetization([0.0,0.0,1.0]) # providing a vector rather than a function

import nfem.visual

nfem.visual.fields2vtkfile([mag.default_simulation_context.field_M],'sphere-initial.vtk',mesh=mag.default_simulation_context.mesh)

mag.advance_time_cvode(0.02)
    
nfem.visual.fields2vtkfile([mag.default_simulation_context.field_M],'sphere-final.vtk',mesh=mag.default_simulation_context.mesh)

