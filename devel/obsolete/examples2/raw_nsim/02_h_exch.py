# Slightly more sophisticated lialg_machine cvode demo script

# XXX NOTE: for now, we also have to pass on dM to the jacobi
# plan. There is a problematic issue with the dof name resolution of
# the EOM LHS! FIXME!

# Notes:
#
# (1) It is somewhat problematic/annoying that the first script does
# not really make clear when a name is a MWE name and when it is a
# field/vector name. Corrected that for this example.

import os,time,sys,math
import nmesh

execfile("../../interface/nsim/linalg_machine.py")

objects=[nmesh.difference(nmesh.ellipsoid([3.0,3.0],[("shift",[0,0])]),
                          [nmesh.ellipsoid([2.5,2.5],[("shift",[0,0])])])
         ]

mesh = nmesh.mesh(objects=objects,
                  a0=0.25,
                  bounding_box=[[-5.0,-5.0],[5.0,5.0]],
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

elem_M=ocaml.make_element("M",[3],2,1)
elem_H=ocaml.make_element("H",[3],2,1)


mwe_M=ocaml.make_mwe("M",raw_mesh,[(1,elem_M)],[])
mwe_H=ocaml.make_mwe("H",raw_mesh,[(1,elem_H)],[])

mwe_dM=ocaml.mwe_sibling(mwe_M,"dM","dM/M",[("M","dM")])
mwe_H_exch=ocaml.mwe_sibling(mwe_H,"H_exch","H_exch/H",[("H","H_exch")])


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
    print "ALPHA ",alpha
    if ix==2:
        return 0.0
    elif ix==0:
        return math.sin(alpha)
    else:
        return -math.cos(alpha)


field_M0=ocaml.raw_make_field(mwe_M,[fun_M0_v2],"","")
field_M=ocaml.raw_make_field(mwe_M,[],"","")

eq_rhs="""
%range i:3, j:3, k:3, p:3, q:3;

dM(i) <- 0.1*eps(i,j,k)*M(j)*H_exch(k) + (-0.3)*eps(i,j,k)*M(j)*eps(k,p,q)*M(p)*H_exch(q);
""";


lam=make_linalg_machine("lam_precession",
                        mwes=[mwe_M,mwe_dM,mwe_H_exch],
                        vectors=[lam_vector(name="v_M",mwe_name="M"),
                                 lam_vector(name="v_dM",mwe_name="dM"),
                                 lam_vector(name="v_H_exch",mwe_name="H_exch"),
                                 ],
                        operators=[lam_operator("op_H_exch","H_exch","M",
                                                "-2.5*<d/dxj H_exch(k)||d/dxj M(k)>, j:2,k:3",
                                                for_jacobian=True),
                                   ],
                        local_operations=[lam_local("set_dM",
                                                    field_mwes=["M","dM","H_exch"],
                                                    equation=eq_rhs),
                                          ],
                        jacobi_plans=[lam_jplan("jacobian", # name
                                                "dM",
                                                ["M","H_exch"], # mwe_names
                                                [[],["op_H_exch"]], # opt_derive_me
                                                eq_rhs   # eom_str
                                                ),
                                      ],
                        programs=[lam_program("set_H_exch",
                                              commands=[["SM*V","op_H_exch","v_M","v_H_exch"],
                                                        ["CFBOX","H_exch","v_H_exch"],
                                                        ["DEBUG","Exchange","v_H_exch",3],
                                                        ]),
                                  lam_program("rhs",
                                              args_fields=[["arg_dM","dM"],["arg_M","M"]], # [arg_name,mwe_name]
                                              commands=[["DISTRIB","arg_M","v_M"],
                                                        ["GOSUB", "set_H_exch"],
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

(cvode,fun_timing)=ocaml.raw_make_linalg_machine_cvode(\
    lam, field_M0,
    "jacobian","execute_jplan",
    "rhs",
    True,2,300
    )

def print_field(field,name="M",pos=[2.8,0.0],t=None):
    probed=ocaml.probe_field(field,name,pos)
    v=probed[0][1]
    if t!=None:
        print "%8.6f  %8.6f %8.6f %8.6f" % (t,v[0],v[1],v[2])
    else:
        print "%8.6f %8.6f %8.6f" % (v[0],v[1],v[2])

print "*** CVODE ***",cvode
sys.stdout.flush()

print_field(field_M0,t=0.0)

for tt in range(1,100):
    target_time=tt*0.1
    ocaml.raw_cvode_advance(cvode,field_M,target_time,-1)
    print_field(field_M,t=target_time)
