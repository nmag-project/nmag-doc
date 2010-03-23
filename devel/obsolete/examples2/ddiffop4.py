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
                             ],
                   a0=0.8,
                   max_steps=400,
                   cache_name="ddiffop4_mesh_rod")

the_mesh = mesh_obj.raw_mesh

print "MESH: ",the_mesh,type(the_mesh)

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
mwe_phi_ddd = ocaml.mwe_sibling(mwe_scalar,"mwe_phi_ddd","phi_ddd/S",[("S","phi_ddd")])
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
    # return 1.0/(1.0+w*w)
    return 1.0

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

field_phi_ddd=ocaml.raw_make_field(mwe_phi_ddd,
                                   [],
                                   "phi_ddd",
                                   "field_phi_ddd")


cofield_rho=ocaml.raw_make_cofield(mwe_rho,
                                   [],
                                   "rho",
                                   "field_rho")

# nfem.field_print_contents(cofield_rho) # initially zero, okay.

#  C_(N|B) = M P_(NB)
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

#print "*** OPERATOR M_NB ***"
#sys.stdout.flush()

op_dirichlet_effective_rho0 = ocaml.fem_operator_from_ddiffop_string(\
    "<d/dxk rho[vol=1]         || d/dxk phi[boundary=1/-2,1/-3]>\
    +<d/dxk rho[boundary=1/-1] || d/dxk phi[boundary=1/-2,1/-3]>;\
    (L||R)=(*||phi[boundary=1/-2,1/-3]),k:3",
    [],mwe_rho,mwe_phi)

#print "*** DONE OPERATOR M_NB ***"
#sys.stdout.flush()

#op_dirichlet_effective_rho = ocaml.fem_operator_from_ddiffop_string(\
#    "<d/dxk rho[vol=1]         |sigma| d/dxk phi[boundary=1/-2,1/-3]>\
#    +<d/dxk rho[boundary=1/-1] |sigma| d/dxk phi[boundary=1/-2,1/-3]>;\
#    (L||R)=(*||phi[boundary=1/-2,1/-3]),k:3",
#    [field_sigma],mwe_rho,mwe_phi)

op_laplace = ocaml.fem_operator_from_ddiffop_string(\
    # "-<d/dxk rho[vol]|sigma|d/dxk phi[vol]>; rho[boundary=1/2,1/3]=phi[boundary=1/2,1/3],k:3",
    "-<d/dxk rho[vol]           |sigma| d/dxk phi[vol]>\
     -<d/dxk rho[boundary=1/-1] |sigma| d/dxk phi[vol]>\
     -<d/dxk rho[vol]           |sigma| d/dxk phi[boundary=1/-1]>\
     -<d/dxk rho[boundary=1/-1] |sigma| d/dxk phi[boundary=1/-1]>\
     ; rho[boundary=1/-2,1/-3]=phi[boundary=1/-2,1/-3],k:3",
    [field_sigma],mwe_rho,mwe_phi)

laplace_solver=ocaml.fem_solver_from_operator(op_laplace)

nfem.field_print_contents(field_phi_dbc)

cofield_rho=ocaml.apply_fem_operator([cofield_rho],
                                     op_dirichlet_effective_rho0,
                                     field_phi_dbc)

#nfem.field_print_contents(cofield_rho)
# XXX The following two lines below do not work - bug in cofield_to_field?!?
field_rho_eff=ocaml.cofield_to_field([],cofield_rho,True)
nfem.field_print_contents(field_rho_eff)

ocaml.apply_fem_solver(cofield_rho,laplace_solver,field_phi)

ocaml.field_push(field_phi,field_phi_ddd) # "save" intermediate field

ocaml.field_push(field_phi_dbc,field_phi) # add proper boundary values

nfem.visual.fields2vtkfile([field_rho_eff,field_phi_ddd,field_phi],'/tmp/ddiffop4.vtk',mesh=mesh_obj)
