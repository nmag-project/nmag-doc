# This example handles both demag and exchange.

import os,time,sys,math
import nmesh

execfile("../../interface/nsim/linalg_machine.py")

#ocaml.init_hlib("/home/fangohr/build/HLib-1.3/Library/.libs/libhmatrix-1.3.so")
ocaml.init_hlib("/home/tf/HLib-1.3/Library/.libs/libhmatrix-1.3.so")

objects=[nmesh.ellipsoid([3.0,3.0,3.0])
         ]


mesh = nmesh.mesh(objects=objects,
                  a0=1.0,
                  bounding_box=[[-5.0,-5.0,-5.0],[5.0,5.0,5.0]],
                  cache_name="geometry.mesh",
                  )

raw_mesh=mesh.raw_mesh

if ocaml.petsc_is_mpi():
    print "*** PARALLEL EXECUTION ***"
    nr_nodes=ocaml.petsc_mpi_nr_nodes()
    nr_points=ocaml.mesh_nr_points(raw_mesh)
    z=nr_points/nr_nodes
    distrib = [int(round(z*(i+1)))-int(round(z*i)) for i in range(0,nr_nodes)]
    slack=nr_points-reduce(lambda x,y:x+y,distrib)
    distrib[0] = distrib[0] + slack
    print "*** RAW MESH %s *** DISTRIB %s ***" %(repr(raw_mesh),repr(distrib))
    ocaml.mesh_set_vertex_distribution(raw_mesh,distrib)

elem_V=ocaml.make_element("V",[3],3,1) # vector
elem_S=ocaml.make_element("S",[],3,1) # scalar element

mwe_V=ocaml.make_mwe("V",raw_mesh,[(1,elem_V)],[]) 
mwe_S=ocaml.make_mwe("S",raw_mesh,[(1,elem_S)],[])


mwe_M=ocaml.mwe_sibling(mwe_V,"M","M/V",[("V","M")])
mwe_dM=ocaml.mwe_sibling(mwe_V,"dM","dM/V",[("V","dM")])
mwe_H_demag=ocaml.mwe_sibling(mwe_V,"H_demag","H_demag/V",[("V","H_demag")])
mwe_H_exch=ocaml.mwe_sibling(mwe_V,"H_exch","H_exch/V",[("V","H_exch")])
mwe_H_total=ocaml.mwe_sibling(mwe_V,"H_total","H_total/V",[("V","H_total")])

mwe_rho=ocaml.mwe_sibling(mwe_S,"rho","rho/S",[("S","rho")])
mwe_phi=ocaml.mwe_sibling(mwe_S,"phi","phi/S",[("S","phi")])

def fun_M0(dof_name,dof_pos):
    ix=dof_name[1][0]
    alpha=math.atan2(dof_pos[0],dof_pos[1])
    if ix==2:
        return 0.0
    elif ix==0:
        return math.sin(alpha)
    else:
        return -math.cos(alpha)

def fun_M0_v2(dof_name,dof_pos):
    ix=dof_name[1][0]
    alpha=math.atan2(dof_pos[0],dof_pos[1])
    alpha=math.exp(-alpha*alpha*3.0)
    alpha=0.0
    print "ALPHA ",alpha
    if ix==2:
        return 0.0
    elif ix==0:
        return math.sin(alpha)
    else:
        return -math.cos(alpha)


field_M0=ocaml.raw_make_field(mwe_M,[fun_M0_v2],"","")
field_M=ocaml.raw_make_field(mwe_M,[],"","")
field_H_demag=ocaml.raw_make_field(mwe_H_demag,[],"","")

eq_rhs="""
%range i:3, j:3, k:3, p:3, q:3;

dM(i) <-
    0.1*eps(i,j,k)*M(j)*H_total(k)
  + (-0.3)*eps(i,j,k)*M(j)*eps(k,p,q)*M(p)*H_total(q)
  + (-1.0)*(M(j)*M(j)+(-1.0))*M(i);
""";

eq_H_total="""
%range j:3;

H_total(j) <- H_demag(j) + H_exch(j);
"""


lam=make_linalg_machine("lam_mumag",
                        mwes=[mwe_M,mwe_dM,mwe_H_demag,mwe_H_exch,mwe_H_total,mwe_rho, mwe_phi],
                        vectors=[lam_vector(name="v_M",mwe_name="M"),
                                 lam_vector(name="v_dM",mwe_name="dM"),
                                 lam_vector(name="v_H_exch",mwe_name="H_exch"),
                                 lam_vector(name="v_H_demag",mwe_name="H_demag"),
                                 lam_vector(name="v_H_total",mwe_name="H_total"),
                                 lam_vector(name="v_rho",mwe_name="rho"),
                                 lam_vector(name="v_rho_s",mwe_name="rho"),
                                 lam_vector(name="v_phi",mwe_name="phi"),
                                 lam_vector(name="v_phi1",mwe_name="phi"),
                                 lam_vector(name="v_phi2",mwe_name="phi"),
                                 lam_vector(name="v_phi1b",mwe_name="phi",restriction="phi[boundary=-1/*]"),
                                 lam_vector(name="v_phi2b",mwe_name="phi",restriction="phi[boundary=-1/*]"),
                                 ],
                        operators=[lam_operator("op_H_exch","H_exch","M",
                                                "-2.5*<d/dxj H_exch(k)||d/dxj M(k)>, j:3,k:3",
                                                for_jacobian=True),
                                   # Note: have to operate on normalized magnetization
                                   # and include magnitude factor! Also have to parametrise properly
                                   # for multiple magnetizations!
                                   lam_operator("op_div_M","rho","M",
                                                "<rho||d/dxj M(j)>+<rho||D/Dxj M(j)>, j:3"),
                                   lam_operator("op_neg_laplace_phi","rho","phi",
                                                "<d/dxj rho || d/dxj phi>;gauge_fix:phi, j:3"
                                                # XXX add matoptions!
                                                ),
                                   lam_operator("op_grad_phi","H_demag","phi",
                                                "-<H_demag(j) || d/dxj phi>, j:3"),
                                   lam_operator("op_laplace_DBC","phi","phi",
                                                "-<d/dxj phi[vol] || d/dxj phi[vol]>;phi[boundary=-1/*]=phi[boundary=-1/*], j:3"),
                                   lam_operator("op_load_DBC","phi","phi",
                                                "<d/dxj phi[vol] || d/dxj phi[boundary=-1/*]>;(L||R)=(*||phi[boundary=-1/*]), j:3"),
                                   ],
                        bem_matrices=[lam_bem(name="BEM",
                                              is_hlib=True,
                                              #is_hlib=False,
                                              mwe_name="phi",
                                              inside_regions=[1],
                                              )],
                        ksps=[lam_ksp(name="solve_neg_laplace_phi",
                                      matrix_name="op_neg_laplace_phi"),
                              lam_ksp(name="solve_laplace_DBC",
                                      matrix_name="op_laplace_DBC"),
                              ],
                        local_operations=[lam_local("addup_H_total",
                                                    field_mwes=["H_total","H_exch","H_demag"],
                                                    equation=eq_H_total),
                                          lam_local("set_dM",
                                                    field_mwes=["M","dM","H_total"],
                                                    equation=eq_rhs),
                                          ],
                        jacobi_plans=[lam_jplan("jacobian", # name
                                                "dM",
                                                ["M","H_total"], # mwe_names
                                                [[],["op_H_exch"]], # opt_derive_me
                                                eq_rhs # eom_str
                                                ),
                                      ],
                        programs=[lam_program("set_H_exch",
                                              commands=[["SM*V","op_H_exch","v_M","v_H_exch"],
                                                        ["CFBOX","H_exch","v_H_exch"],
                                                        ["DEBUG","Exchange","v_H_exch",3],
                                                        ]),
                                  lam_program("set_H_demag",
                                              commands=[["SM*V","op_div_M","v_M","v_rho"],
                                                        ["SCALE","v_rho",-1.0],
                                                        ["SOLVE","solve_neg_laplace_phi","v_rho","v_phi1"],
                                                        ["PULL-FEM","phi","phi[boundary=-1/*]","v_phi1","v_phi1b"],
                                                        ["DM*V","BEM","v_phi1b","v_phi2b"],
                                                        ["SM*V","op_load_DBC","v_phi2b","v_rho_s"],
                                                        ["SOLVE","solve_laplace_DBC","v_rho_s","v_phi2"],
                                                        ["PUSH-FEM","phi","phi[boundary=-1/*]","v_phi2b","v_phi2"],
                                                        ["AXPBY",1.0,"v_phi1",0.0,"v_phi"],
                                                        ["AXPBY",1.0,"v_phi2",1.0,"v_phi"],
                                                        ["SM*V","op_grad_phi","v_phi","v_H_demag"],
                                                        ["CFBOX","H_demag","v_H_demag"]
                                                        ]
                                              ),
                                  lam_program("set_H_total",
                                              commands=[
                                                         ["SITE-WISE","addup_H_total",["v_H_total","v_H_exch","v_H_demag"],[],[]]
                                                         ]),
                                  lam_program("rhs",
                                              args_fields=[["arg_dM","dM"],["arg_M","M"]], # [arg_name,mwe_name]
                                              commands=[["DISTRIB","arg_M","v_M"],
                                                        ["GOSUB", "set_H_exch"],
                                                        ["GOSUB", "set_H_demag"],
                                                        ["GOSUB", "set_H_total"],
                                                        ["SITE-WISE","set_dM",["v_M","v_dM","v_H_exch"],[],[]],
                                                        ["COLLECT","arg_dM","v_dM"]
                                                        ]),
                                  lam_program("execute_jplan",
                                              commands=[["GOSUB", "set_H_exch"],
                                                        ["JPLAN","jacobian",["v_M","v_H_exch"]]
                                                        ]),
                                  ]
                        )

print "*** LAM ***",lam
sys.stdout.flush()

cvode=ocaml.raw_make_linalg_machine_cvode(\
    lam, field_M0,
    "jacobian","execute_jplan",
    "rhs"
    )

def print_field(field,name="M",pos=[1.0,0.0,0.0],t=None):
    probed=ocaml.probe_field(field,name,pos)
    v=probed[0][1]
    if t!=None:
        print "%8.6f  %8.6f %8.6f %8.6f" % (t,v[0],v[1],v[2])
    else:
        print "%8.6f %8.6f %8.6f" % (v[0],v[1],v[2])

print "*** CVODE ***",cvode
sys.stdout.flush()

print_field(field_M0,t=0.0)


import nfem.visual

for tt in range(1,2):
    target_time=tt*0.1
  
    ocaml.raw_cvode_advance(cvode,field_M,target_time,-1)
    
    ocaml.linalg_machine_get_field(lam,field_H_demag,"v_H_demag")

    nfem.visual.fields2vtkfile([field_M,field_H_demag],'test.vtk', mesh)
    print_field(field_M,t=target_time)




