#!../bin/nmesh2

import nmesh, Numeric, math,time, sys
import nfem
import nfem.visual

nfem.set_default_dimension(2)
nfem.set_default_order(1)

# Simulation parameters
sigma0=1.0

print
"""
** Example: meshing half ring and compute current density for contacts at either side.
This is for possible collaboration with Uni Hamburg. **
"""

##### Creating the mesh #####

ring = nmesh.difference(nmesh.ellipsoid([4.0,4.0]),[nmesh.ellipsoid([2.5,2.5])])

halfring=nmesh.intersect([ring,nmesh.box([-4.0,-0.0],[4.0,4.0])])

the_mesh = nmesh.mesh(objects = [halfring],
                      cache_name="halfring",
                      a0=0.2,
                      bounding_box=[[-4.0,-4],[4.0,4.0]],
                      neigh_force_scale = 1.,
                      initial_settling_steps = 50,
                      max_relaxation = 4,
                      max_steps=400
                      )

nfem.set_default_mesh(the_mesh)

the_mesh.save("themesh.nmesh")

##### Making the elements... #####

# conductivity (scalar)
element_sigma     = nfem.make_element("sigma",[]);
element_drho_by_dt= nfem.make_element("drho_by_dt",[]);
element_phi       = nfem.make_element("phi",[]);
element_J         = nfem.make_element("J",[2]);

mwe_sigma      = nfem.make_mwe("mwe_sigma",     [(1,element_sigma)])
mwe_drho_by_dt = nfem.make_mwe("mwe_drho_by_dt",[(1,element_drho_by_dt)])
mwe_phi        = nfem.make_mwe("mwe_phi",       [(1,element_phi)])
mwe_J          = nfem.make_mwe("mwe_J",         [(1,element_J)])

diffop_laplace=nfem.diffop("-<d/dxj drho_by_dt|sigma|d/dxj phi>, j:2")
diffop_J=nfem.diffop("<J(k)|sigma|d/dxk phi>, k:2")

# Initial conductivity is spatially constant:

# Later on, we will modify this field:



def fun_sigma0(dof_name_indices,position):
    return sigma0

field_sigma=nfem.make_field(mwe_sigma,fun_sigma0)

def fun_sigma0(dof_name_indices,position):
    return sigma0

# Later on, we will modify this field:

field_sigma=nfem.make_field(mwe_sigma,fun_sigma0)

# Dirichlet Boundary Conditions on our sample:

def laplace_dbc_values(dof_name_indices,coords):
    x,y=coords
    if x > 0:
        if y < 0.0+0.05:
            return 0.0
    if x < 0:
        if y < 0.0+0.05:
            return -1.0
    else:
        print "unused laplace_dbc at x =",x,"and y =",y


#this line causes a crash 
#def laplace_dbc(dof_name_indices,coord):

##where as this one works.
def laplace_dbc(coords):

    #the problem is the wrong signature of the function.
    #Can we catch this?
    print coords
    x,y=coords
    if x > 0:
        if y < 0.0+0.05:
            print "Allocate  dbc at (%f,%f)" % (x,y)
            return 1
    elif x < 0:
        if y < 0.0+0.05:
            print "Allocate  dbc at (%f,%f)" % (x,y)
            return 1
    else:
        return 0


       

cofield_drho_by_dt=nfem.make_cofield(mwe_drho_by_dt)

prematrix_laplace=nfem.prematrix(diffop_laplace,
                                  mwe_drho_by_dt,mwe_phi,
                                  mwe_mid=mwe_sigma)

prematrix_J=nfem.prematrix(diffop_J,
                                   mwe_J,mwe_phi,
                                   mwe_mid=mwe_sigma)

compute_J=nfem.prematrix_applicator(prematrix_J,
                                     field_mid=field_sigma)
print "4"
laplace_solver=nfem.laplace_solver(prematrix_laplace,
                                    dirichlet_bcs=[(-1,1,laplace_dbc)],
                                    mwe_mid=field_sigma)

field_phi = laplace_solver(cofield_drho_by_dt,
                           dbc_values=laplace_dbc_values)

field_J = nfem.cofield_to_field(compute_J(field_phi))

nfem.visual.fields2vtkfile([field_phi,field_J],'result.vtk',the_mesh)

print "Try to visualise with "
print "mayavi -d run_fem-halfring/result.vtk -m SurfaceMap -m VelocityVector"

print "The current density magnitude can be calculated in MayaVi (Filter: extractVectorNorm)"


