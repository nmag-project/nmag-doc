#!../bin/nmesh2


import math,time, sys,logging
import nfem,nmesh

ocaml.sys_profiling("on")

#get user-specific logger
log=logging.getLogger('user')
#could modify debug level for this logger
log.setLevel(logging.DEBUG)


nfem.set_default_dimension(3)
nfem.set_default_order(1)

##### Creating the mesh #####

# For now, we use a very very simple mesh...

ball = nmesh.ellipsoid([3.0,3.0,3.0])
density = "density=1.;"


the_mesh = nmesh.mesh(objects = [ball],
                      cache_name="bem-ball",
                      a0=1.0,
                      bounding_box=[[-5.0,-5.0,-5.0],[5.0,5.0,5.0]],
                      neigh_force_scale = 1.,
                      density = density,
                      initial_settling_steps = 50,
                      max_relaxation = 4,
                      max_steps=500
                      )

nfem.set_default_mesh(the_mesh)

##### Making the elements... #####

element_M     = nfem.make_element("M",[3]);
element_H     = nfem.make_element("H",[3]);
element_rho_M = nfem.make_element("rho_M",[]);
element_phi_M = nfem.make_element("phi_M",[]);

mwe_M      = nfem.make_mwe("mwe_M",     [(1,element_M)])
mwe_H      = nfem.make_mwe("mwe_H",     [(1,element_H)])
mwe_rho_M  = nfem.make_mwe("mwe_rho_M",     [(1,element_rho_M)])
mwe_phi_M  = nfem.make_mwe("mwe_phi_M",     [(1,element_phi_M)])

log.info("OK 1!")

# Initial magnetization is spatially constant: M=(1 0 0)

#def fun_M0(dof_name_indices,pos):
#    r=math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])
#    if dof_name_indices[1][0] == 0 and r < 2.9:
#        return 1.0
#    else:
#        return 0.0

def fun_M0(dof_name_indices,pos):
    if dof_name_indices[1][0] == 0:
        return 1.0
    else:
        return 0.0


diffop_laplace=nfem.diffop("-<d/dxj rho_M||d/dxj phi_M>, j:3")
diffop_div_M=nfem.diffop("<rho_M||d/dxj M(j)>, j:3")
diffop_grad_phi=nfem.diffop("<H(j)||d/dxj phi_M>, j:3")

log.info("OK 2!")

#This is " $\ laplace -phi = =rho_M$
prematrix_laplace=nfem.prematrix(diffop_laplace,mwe_rho_M,mwe_phi_M)
prematrix_div_M=nfem.prematrix(diffop_div_M,mwe_rho_M,mwe_M,ignore_jumps=False)
prematrix_grad_phi=nfem.prematrix(diffop_grad_phi,mwe_H,mwe_phi_M)

log.info("OK 3!")

solve_bem=nfem.laplace_solver_bem(prematrix_laplace,inside_regions=[1])

log.info("OK 4!")

compute_div_M=nfem.prematrix_applicator(prematrix_div_M,interface_coeffs=[(1,-1,1.0)])
compute_grad_phi=nfem.prematrix_applicator(prematrix_grad_phi)

log.info("OK 5!")

field_M0=nfem.make_field(mwe_M,fun_M0)

cofield_div_M=compute_div_M(field_M0)

# DDD
ddd_field_div_M=nfem.cofield_to_field(cofield_div_M)
print "DDD div M: ",nfem.probe_field(ddd_field_div_M,[2.7,0.2,0.0]) # we may add: ,name_stem="rho_M")
# END DDD

field_phi_M=solve_bem(cofield_div_M)

log.info("OK 6!")
field_H=nfem.cofield_to_field(compute_grad_phi(field_phi_M))

log.info("OK 7!")
# Also, we do direct summation, just for comparison:

direct_summer = ocaml.dipole_field_direct_summer(1, # nr cpus (processes to fork to)
                                                 0.0, # minimal distance
                                                 mwe_H,
                                                 mwe_M,
                                                 "H","M")

field_H_dsum=nfem.make_field(mwe_H)

direct_summer(field_H_dsum,field_M0)

for i in range(0,10):
    pos=[-1.5+i*0.2,0.0,0.0]
    print "Position ",pos," H: ",nfem.probe_field(field_H,pos), "H_dsum: ",nfem.probe_field(field_H_dsum,pos)

print "*** PROFILING ***\n",ocaml.sys_profiling("report")


sys.exit(0)

#create vtk files with demag field and potential
import nfem.visual as v
#v.scalarfield2vtkfile(field_phi_M,"phi_M",'demagphi.vtk',the_mesh)
#v.vectorfield2vtkfile(field_H,"H",'demagH.vtk',the_mesh)
#print "Data can be visualized with this command:"
#print "'mayavi -d   demagH.vtk -m VelocityVector -d demagphi.vtk -m SurfaceMap'"
#
v.fields2vtkfile( field_phi_M, 'test.vtk',path='.' )

