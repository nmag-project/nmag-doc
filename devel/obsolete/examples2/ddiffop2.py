# Test example 2: heat flow with spatially changing heat conductivity
#
# Using the discretized differential operator machinery
# from within python...
#

import nmag2 as nmag, nmesh as nm, sys, math, time
import nfem
import nfem.visual

# NOTE: CSG semantics sucks here!

# temporarily commented out by turning this into a string...
"""
the_mesh = nm.mesh(([-5.0,-5.0,-5.0],[5.0,5.0,5.0]), # bounding box
                   objects= \
                    [nm.difference( \
                      nm.intersect([nm.conic([0.0,0.0,-2.0],7.0,
                                             [0.0,0.0,2.0],7.0),
                                    nm.box([-4.0,-4.0,-4.0],[4.0,4.0,4.0])]))],
                   a0=1.2,
                   max_steps=500,
                   cache_name="ddiffop_mesh")
"""

mesh_obj = nm.mesh(([-5.0,-5.0,-5.0],[5.0,5.0,5.0]), # bounding box
                   objects= [nm.conic([0.0,0.0,-2.0],4.0,[0.0,0.0,2.0],4.0)],
                   a0=1.2,
                   max_steps=500,
                   cache_name="ddiffop_mesh")

the_mesh = mesh_obj.raw_mesh

print "MESH: ",the_mesh,type(the_mesh)


elem_T = ocaml.make_element("T",[],3,1)
elem_sigma = ocaml.make_element("sigma",[],3,1)
elem_j_q = ocaml.make_element("j_q",[3],3,1)

mwe_T=ocaml.make_mwe("mwe_T",the_mesh,[(1,elem_T)],[])
mwe_sigma=ocaml.make_mwe("mwe_sigma",the_mesh,[(1,elem_sigma)],[])
mwe_j_q=ocaml.make_mwe("mwe_j_q",the_mesh,[(1,elem_j_q)],[])

def initial_T(dof,coords):
    # temp = math.sin(2.0*coords[1])
    temp = coords[1]
    return temp

def initial_sigma(dof,coords):
    # temp = math.sin(2.0*coords[1])
    return 1.0/(1.0+0.9*math.sin(2.0*coords[0]))


field_T=ocaml.raw_make_field(mwe_T,
                             [initial_T],
                             "",
                             "field_T")

field_sigma=ocaml.raw_make_field(mwe_sigma,
                                 [initial_sigma],
                                 "",
                                 "field_sigma")

field_j_q=ocaml.raw_make_field(mwe_j_q,
                               [],
                               "",
                               "field_j_q")

op_grad_T = ocaml.fem_operator_from_ddiffop_string("<j_q(k)|sigma|d/dxk T>,k:3",[field_sigma],mwe_j_q,mwe_T)

field_j_q=ocaml.apply_fem_operator_ff([field_j_q],op_grad_T,field_T)

nfem.visual.fields2vtkfile([field_j_q,field_T],'/tmp/ddiffop2.vtk',mesh=mesh_obj)
