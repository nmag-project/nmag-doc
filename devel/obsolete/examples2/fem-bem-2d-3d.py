#!../bin/nmesh2


import math,time, sys,logging
import nfem,nmesh

#get user-specific logger
log=logging.getLogger('user')
#could modify debug level for this logger
log.setLevel(logging.DEBUG)


nfem.set_default_dimension(2)
nfem.set_default_order(1)

thickness2d=0.4

# Initial magnetization is spatially constant:
# (We use the same M0 for both cases!)

def fun_M0(dof_name_indices,pos):
    m=[0.0,1.0,0.0]
    return m[dof_name_indices[1][0]]

# Note that we can easily also use the 3d differential operators for
# the 2d problems...

diffop_laplace=nfem.diffop("-<d/dxj rho_M||d/dxj phi_M>, j:3")
diffop_div_M=nfem.diffop("<rho_M||d/dxj M(j)>, j:3")
diffop_grad_phi=nfem.diffop("<H(j)||d/dxj phi_M>, j:3")


##### Creating the mesh #####

# For now, we use a very very simple mesh...

double_disc_2d = nmesh.union([nmesh.ellipsoid([3.0,3.0], transform=[("shift",[-1.0,0.0])]),
                              nmesh.ellipsoid([3.0,3.0], transform=[("shift",[1.0,0.0])])])

double_disc_3d = nmesh.union([nmesh.conic([-1.0,0.0,thickness2d*0.5],3.0,[-1.0,0.0,-thickness2d*0.5],3.0),
                              nmesh.conic([ 1.0,0.0,thickness2d*0.5],3.0,[ 1.0,0.0,-thickness2d*0.5],3.0)])

density = "density=1.;"


mesh_2d = nmesh.mesh(objects = [double_disc_2d],
                     cache_name="double-disc-2d",
                     a0=0.4,
                     bounding_box=[[-5.0,-5.0],[5.0,5.0]],
                     neigh_force_scale = 1.,
                     density = density,
                     initial_settling_steps = 50,
                     max_relaxation = 4,
                     max_steps=500
                     )

mesh_3d = nmesh.mesh(objects = [double_disc_3d],
                     cache_name="double-disc-3d",
                     a0=0.4,
                     bounding_box=[[-5.0,-5.0,-0.3],[5.0,5.0,0.3]],
                     neigh_force_scale = 1.,
                     density = density,
                     initial_settling_steps = 50,
                     max_relaxation = 4,
                     max_steps=500
                     )

# === The 2.5d computation ===

nfem.set_default_mesh(mesh_2d)

##### Making the elements... #####

elem2d_M     = nfem.make_element("M",[3],dim=2);
elem2d_H     = nfem.make_element("H",[3],dim=2);
elem2d_rho_M = nfem.make_element("rho_M",[],dim=2);
elem2d_phi_M = nfem.make_element("phi_M",[],dim=2);

mwe2d_M      = nfem.make_mwe("mwe2d_M",     [(1,elem2d_M)])
mwe2d_H      = nfem.make_mwe("mwe2d_H",     [(1,elem2d_H)])
mwe2d_rho_M  = nfem.make_mwe("mwe2d_rho_M", [(1,elem2d_rho_M)])
mwe2d_phi_M  = nfem.make_mwe("mwe2d_phi_M", [(1,elem2d_phi_M)])

#This is " $\ laplace -phi = =rho_M$
pmx2d_laplace=nfem.prematrix(diffop_laplace,mwe2d_rho_M,mwe2d_phi_M)
pmx2d_div_M=nfem.prematrix(diffop_div_M,mwe2d_rho_M,mwe2d_M,ignore_jumps=False)
pmx2d_grad_phi=nfem.prematrix(diffop_grad_phi,mwe2d_H,mwe2d_phi_M)

solve_bem_2d=nfem.laplace_solver_bem(pmx2d_laplace,inside_regions=[1], thickness=thickness2d)

compute_div_M_2d=nfem.prematrix_applicator(pmx2d_div_M,interface_coeffs=[(1,-1,1.0)])
compute_grad_phi_2d=nfem.prematrix_applicator(pmx2d_grad_phi)

field2d_M0=nfem.make_field(mwe2d_M,fun_M0)

cofield2d_div_M=compute_div_M_2d(field2d_M0)

field2d_phi_M=solve_bem_2d(cofield2d_div_M)

field2d_H=nfem.cofield_to_field(compute_grad_phi_2d(field2d_phi_M))

# DDDDDD

print "phi_M left: ",nfem.probe_field(field2d_phi_M,[-3.8,0.0])
print "phi_M right: ",nfem.probe_field(field2d_phi_M,[3.8,0.0])

nfem.plot_scalar_field(field2d_phi_M,
                        "phi_M","/tmp/plot-phi-M.ps",
                        plot_edges=False,
                        color_scheme=[(-0.3,[0.1,0.1,1.0]),
                                      (0.0,[0.1,1.0,0.1]),
                                      (0.3,[1.0,0.1,0.1])])

nfem.plot_scalar_field(nfem.cofield_to_field(cofield2d_div_M),
                        "rho_M","/tmp/plot-rho-M.ps",
                        plot_edges=False,
                        color_scheme=[(-2.0,[0.1,0.1,1.0]),
                                      (0.0,[0.1,1.0,0.1]),
                                      (2.0,[1.0,0.1,0.1])])


# End DDDDDD

# we still have to correct the Z-field component. For now, we do this
# with a site-wise applicator...

adjust_H_z_2d=nfem.site_wise_applicator(parameter_names=[],
                                         # all the names of extra parameters
                                         code="H(2)=-M(2);",
                                         field_mwes=[mwe2d_M,mwe2d_H]
                                         )
adjust_H_z_2d([],fields=[field2d_M0,field2d_H])

# nfem.field_print_contents(field2d_H) # DDD

print "*** The 2.5-dimensional result ***"

for i in range(0,76):
    pos=[-3.8+i*0.1,0.0]
    print "Position ",pos," H: ",nfem.probe_field(field2d_H,pos)

# === The 3d computation ===

nfem.set_default_mesh(mesh_3d)

##### Making the elements... #####

elem3d_M     = nfem.make_element("M",[3],dim=3);
elem3d_H     = nfem.make_element("H",[3],dim=3);
elem3d_rho_M = nfem.make_element("rho_M",[],dim=3);
elem3d_phi_M = nfem.make_element("phi_M",[],dim=3);

mwe3d_M      = nfem.make_mwe("mwe3d_M",     [(1,elem3d_M)])
mwe3d_H      = nfem.make_mwe("mwe3d_H",     [(1,elem3d_H)])
mwe3d_rho_M  = nfem.make_mwe("mwe3d_rho_M",     [(1,elem3d_rho_M)])
mwe3d_phi_M  = nfem.make_mwe("mwe3d_phi_M",     [(1,elem3d_phi_M)])

pmx3d_laplace=nfem.prematrix(diffop_laplace,mwe3d_rho_M,mwe3d_phi_M)
pmx3d_div_M=nfem.prematrix(diffop_div_M,mwe3d_rho_M,mwe3d_M,ignore_jumps=False)
pmx3d_grad_phi=nfem.prematrix(diffop_grad_phi,mwe3d_H,mwe3d_phi_M)

solve_bem_3d=nfem.laplace_solver_bem(pmx3d_laplace,inside_regions=[1])

compute_div_M_3d=nfem.prematrix_applicator(pmx3d_div_M,interface_coeffs=[(1,-1,1.0)])
compute_grad_phi_3d=nfem.prematrix_applicator(pmx3d_grad_phi)

field3d_M0=nfem.make_field(mwe3d_M,fun_M0)

cofield3d_div_M=compute_div_M_3d(field3d_M0)

field3d_phi_M=solve_bem_3d(cofield3d_div_M)

field3d_H=nfem.cofield_to_field(compute_grad_phi_3d(field3d_phi_M))

print "*** The 3-dimensional result ***"

for i in range(0,76):
    pos=[-3.8+i*0.1,0.0,0.0]
    print "Position ",pos," H: ",nfem.probe_field(field3d_H,pos)
