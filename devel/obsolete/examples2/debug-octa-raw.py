#
# (C) 2006 Dr. Thomas Fischbacher
# Checking bulk and surface contribs for a homogeneously magnetized octahedron.


import nfem
import sys,math,time
import ocaml #need this as long as we use ocaml.probe_field

mesh = ocaml.mesh_readfile("./debug-octa.mesh")

print "MESH: ",mesh

elem_m=ocaml.make_element("m_X",[3],3,1)
elem_H_demag=ocaml.make_element("H_demag",[3],3,1)

mwe_m=ocaml.make_mwe("mwe_m",mesh,[(1,elem_m)])
mwe_h=ocaml.make_mwe("mwe_h",mesh,[(1,elem_H_demag)])

def initial_m(dof,pos): # radially outward, zero at center
    ix=dof[1][0]
    return pos[ix]

field_m=ocaml.raw_make_field(mwe_m,[initial_m],"")

make_field_h=ocaml.ddd_demag_fun_3d("<S||d/dxj m_X(j)>, j:3",mwe_h,mwe_m)

field_h=make_field_h(field_m)

print field_h

sys.exit()

