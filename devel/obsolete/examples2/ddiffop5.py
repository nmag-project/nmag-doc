# Test example 2: heat flow with spatially changing heat conductivity
#
# Using the discretized differential operator machinery
# from within python...
#

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual


mesh_obj = nm.mesh(([-6.0,-6.0,-6.0],[6.0,6.0,6.0]), # bounding box
                   objects= [nm.conic([0.0,0.0,-4.0],2.0,[0.0,0.0,4.0],2.0),
                             ],
                   a0=0.7,
                   max_steps=500,
                   cache_name="ddiffop5_mesh_rod")

the_mesh = mesh_obj.raw_mesh

elem_scalar = ocaml.make_element("S",[],3,1)
elem_vector = ocaml.make_element("V",[3],3,1)

def surface_type(pos):
    if pos[2]<-3.99:
        return -2
    elif pos[2]>3.99:
        return -3
    else:
        return -1

mwe_scalar=ocaml.make_mwe("mwe_S",the_mesh,[(1,elem_scalar),(2,elem_scalar),(3,elem_scalar)],[surface_type])
mwe_vector=ocaml.make_mwe("mwe_V",the_mesh,[(1,elem_vector),(2,elem_vector),(3,elem_vector)],[surface_type])

mwe_phi = ocaml.mwe_sibling(mwe_scalar,"mwe_phi","phi/S",[("S","phi")])
mwe_rho_ddd = ocaml.mwe_sibling(mwe_scalar,"mwe_rho_ddd","rho_ddd/S",[("S","rho_ddd")])
mwe_rho = ocaml.mwe_sibling(mwe_scalar,"mwe_rho","rho/S",[("S","rho")])
mwe_sigma = ocaml.mwe_sibling(mwe_scalar,"mwe_sigma","sigma/S",[("S","sigma")])

mwe_j=ocaml.mwe_sibling(mwe_vector,"mwe_j","j/V",[("V","j")])

def initial_phi_top_bottom(dof,coords):
    if(coords[2]>0):
        return 1.0
    else:
        return -1.0

def initial_sigma(dof,coords):
    x=coords[0]-3.0
    y=coords[1]
    z=coords[2]
    w = x*x+y*y+z*z
    sigma=1.0
    if w<9.0:
        sigma=20.0
    return sigma


field_sigma=ocaml.raw_make_field(mwe_sigma,
                                 [initial_sigma],
                                 "",
                                 "field_sigma")

field_phi_dbc=ocaml.raw_make_field(mwe_phi,
                                   [initial_phi_top_bottom],
                                   "phi[boundary=1/-2,1/-3]",
                                   "field_phi_dbc")

field_phi=ocaml.raw_make_field(mwe_phi,
                               [],
                               "phi",
                               "field_phi")

field_rho_ddd=ocaml.raw_make_field(mwe_rho_ddd,
                                   [],
                                   "rho_ddd",
                                   "field_rho_ddd")


cofield_rho=ocaml.raw_make_cofield(mwe_rho,
                                   [],
                                   "rho",
                                   "field_rho")

# nfem.field_print_contents(cofield_rho) # initially zero, okay.

#  C_(N|B) = M P_(N|B) (Charge_(Nonboundary|Boundary) = Matrix * Potential_(N|B)
#
#  C_N = M_NB P_B + M_NN P_N
#  C_B = M_BB P_B + M_BN P_N
#
#  Known: P_B, C_N
#
#  Want: P_N, C_B
#
#  Hence,
#
#  M_NN P_N = C_N - M_NB P_B

# This is -M_NB P_B:
op_dirichlet_effective_rho = ocaml.fem_operator_from_ddiffop_string(\
    "<d/dxk rho[vol=1]         |sigma| d/dxk phi[boundary=1/-2,1/-3]>\
    +<d/dxk rho[boundary=1/-1] |sigma| d/dxk phi[boundary=1/-2,1/-3]>;\
    (L||R)=(*||phi[boundary=1/-2,1/-3]),k:3",
    [field_sigma],mwe_rho,mwe_phi)

op_laplace = ocaml.fem_operator_from_ddiffop_string(\
    "-<d/dxk rho[vol]           |sigma| d/dxk phi[vol]>\
     -<d/dxk rho[boundary=1/-1] |sigma| d/dxk phi[vol]>\
     -<d/dxk rho[vol]           |sigma| d/dxk phi[boundary=1/-1]>\
     -<d/dxk rho[boundary=1/-1] |sigma| d/dxk phi[boundary=1/-1]>\
     ; rho[boundary=1/-2,1/-3]=phi[boundary=1/-2,1/-3],k:3",
    [field_sigma],mwe_rho,mwe_phi)

laplace_solver=ocaml.fem_solver_from_operator(op_laplace)

cofield_rho=ocaml.apply_fem_operator([cofield_rho],
                                     op_dirichlet_effective_rho,
                                     field_phi_dbc)

ocaml.apply_fem_solver(cofield_rho,laplace_solver,field_phi)
ocaml.field_push(field_phi_dbc,field_phi) # add proper boundary values

field_rho=ocaml.cofield_to_field([],cofield_rho,True)

nfem.visual.fields2vtkfile([field_phi,field_rho,field_sigma],'/tmp/ddiffop5b.vtk',mesh=mesh_obj)
