# Test example 2: heat flow with spatially changing heat conductivity
#
# Using the discretized differential operator machinery
# from within python...
#

# NOTE: Mesher (version from 30-04-2007) has difficulty here, cannot
# generate the mesh. Reverting to the mesher version from 01-04-2007 works.

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual


mesh_obj = nm.mesh(([-6.0,-6.0,-6.0],[6.0,6.0,6.0]), # bounding box
                   objects= [nm.conic([0.0,0.0,-4.0],2.0,[0.0,0.0,4.0],2.0),
                             nm.conic([0.0,0.0,-5.0],2.0,[0.0,0.0,-4.0],2.0),
                             nm.conic([0.0,0.0,5.0],2.0,[0.0,0.0,4.0],2.0)
                             ],
                   a0=0.6,
                   max_steps=100,
                   cache_name="ddiffop_mesh_rod")

the_mesh = mesh_obj.raw_mesh

print "MESH: ",the_mesh,type(the_mesh)


elem_scalar = ocaml.make_element("S",[],3,1)
elem_vector = ocaml.make_element("V",[3],3,1)

mwe_scalar=ocaml.make_mwe("mwe_S",the_mesh,[(1,elem_scalar),(2,elem_scalar),(3,elem_scalar)],[])
mwe_vector=ocaml.make_mwe("mwe_V",the_mesh,[(1,elem_vector),(2,elem_vector),(3,elem_vector)],[])

mwe_phi = ocaml.mwe_sibling(mwe_scalar,"mwe_phi","phi/S",[("S","phi")])
mwe_rho = ocaml.mwe_sibling(mwe_scalar,"mwe_rho","rho/S",[("S","rho")])
mwe_sigma = ocaml.mwe_sibling(mwe_scalar,"mwe_sigma","sigma/S",[("S","sigma")])

mwe_j=ocaml.mwe_sibling(mwe_vector,"mwe_j","j/V",[("V","j")])

def initial_phi_top_bottom(dof,coords):
    if(coords[2]>0):
        return 1.0
    else:
        return -1.0

def initial_sigma(dof,coords):
    w = coords[2]*coords[1]
    return 1.0/(1.0+w*w)

field_sigma=ocaml.raw_make_field(mwe_sigma,
                                 [initial_sigma],
                                 "",
                                 "field_sigma")

field_phi_dbc=ocaml.raw_make_field(mwe_phi,
                                   [initial_phi_top_bottom],
                                   "phi[boundary=1/2,1/3]",
                                   "field_phi_dbc")

field_phi=ocaml.raw_make_field(mwe_phi,
                               [],
                               "phi",
                               "field_phi")

cofield_rho=ocaml.raw_make_cofield(mwe_rho,
                                   [],
                                   "rho",
                                   "field_rho")


op_dirichlet_rho_boundary = ocaml.fem_operator_from_ddiffop_string(\
    "<d/dxk rho|sigma|d/dxk phi[boundary=1/2,1/3]>;(L||R)=(*||phi[boundary=1/2,1/3]),k:3",
    [field_sigma],mwe_rho,mwe_phi)

op_laplace = ocaml.fem_operator_from_ddiffop_string(\
    # "-<d/dxk rho[vol]|sigma|d/dxk phi[vol]>; rho[boundary=1/2,1/3]=phi[boundary=1/2,1/3],k:3",
    "-<d/dxk rho[vol]|sigma|d/dxk phi[vol]>; rho[boundary]=phi[boundary],k:3",
    [field_sigma],mwe_rho,mwe_phi)

laplace_solver=ocaml.fem_solver_from_operator(op_laplace)

cofield_rho=ocaml.apply_fem_operator([cofield_rho],
                                     op_dirichlet_rho_boundary,
                                     field_phi_dbc)

ocaml.apply_fem_solver(cofield_rho,laplace_solver,field_phi)

nfem.visual.fields2vtkfile([field_phi],'/tmp/ddiffop3.vtk',mesh=mesh_obj)
