# This example demonstrates periodic exchange.

import os,time,sys,math
import nmesh
import ocaml
import nfem, nfem.visual

execfile("../../interface/nsim/linalg_machine.py")

#ocaml.init_hlib("/home/fangohr/build/HLib-1.3/Library/.libs/libhmatrix-1.3.so")
ocaml.init_hlib("/home/tf/HLib-1.3/Library/.libs/libhmatrix-1.3.so")

#raw_mesh=ocaml.mesh_readfile("periodic/periodic.nmesh",False)
mesh=nmesh.load_ascii("periodic/periodic-1d-2.nmesh")
raw_mesh=mesh.raw_mesh

print "MESH: ",raw_mesh
sys.stdout.flush()

if ocaml.petsc_is_mpi():
    print "*** PARALLEL EXECUTION *** (not yet supported here!)"
    sys.exit(0) # DDD add parallel support!

elem_phi=ocaml.make_element("phi",[],1,1) # scalar element

mwe_phi=ocaml.make_mwe("phi",raw_mesh,[(0,ocaml.empty_element),(1,elem_phi)],[])
mwe_dphi_dt=ocaml.mwe_sibling(mwe_phi,"dphi_dt","dphi_dt/phi",[("phi","dphi_dt")])
mwe_laplace_phi=ocaml.mwe_sibling(mwe_phi,"laplace_phi","laplace_phi/phi",[("phi","laplace_phi")])


def fun_phi0(dof_name,dof_pos):
    x=dof_pos[0]/10.0
    # sx = math.sin((x-0.4)*2*3.1415926535)
    sx = math.sin((x-0.4)*2*3.1415926535)
    phi=math.exp(-sx*sx)
    return phi

field_phi0=ocaml.raw_make_field(mwe_phi,[fun_phi0],"","")
field_phi=ocaml.raw_make_field(mwe_phi,[],"","")

eq_rhs="dphi_dt <- laplace_phi;"

lam=make_linalg_machine("lam_phi",
                        mwes=[mwe_phi,mwe_dphi_dt,mwe_laplace_phi],
                        vectors=[lam_vector(name="v_phi",mwe_name="phi"),
                                 lam_vector(name="v_dphi_dt",mwe_name="dphi_dt"),
                                 lam_vector(name="v_laplace_phi",mwe_name="laplace_phi"),
                                 ],
                        operators=[lam_operator("op_laplace","laplace_phi","phi",
                                                "-<d/dx0 laplace_phi|| d/dx0 phi>; periodic:laplace_phi",
                                                for_jacobian=True),
                                   ],
                        local_operations=[lam_local("set_dphi_dt",
                                                    field_mwes=["dphi_dt","laplace_phi"],
                                                    equation=eq_rhs),
                                          ],
                        jacobi_plans=[lam_jplan("jacobian", # name
                                                "dphi_dt",
                                                ["phi","laplace_phi"], # mwe_names
                                                [[],["op_laplace"]], # opt_derive_me
                                                eq_rhs # eom_str
                                                ),
                                      ],
                        programs=[lam_program("rhs",
                                              args_fields=[["arg_dphi_dt","dphi_dt"],["arg_phi","phi"]],
                                              commands=[
                                                        ["DISTRIB","arg_phi","v_phi"],
                                                        ["SM*V","op_laplace","v_phi","v_laplace_phi"],
                                                        ["CFBOX","laplace_phi","v_laplace_phi"],
                                                        ["SITE-WISE","set_dphi_dt",["v_dphi_dt","v_laplace_phi"],[],[]],
                                                        ["DEBUG","dphi/dt","v_dphi_dt",11],
                                                        ["COLLECT","arg_dphi_dt","v_dphi_dt"]
                                                        ]),
                                  lam_program("execute_jplan",
                                              commands=[# Note that this Jacobian actually does not
                                                        # depend on laplace_phi,
                                                        # so we need not make sure it is up-to-date.
                                                        # Nevertheless, we have to provide it.
                                                         ["JPLAN","jacobian",["v_phi","v_laplace_phi"]]
                                                        ]),
                                  ]
                        )

print "*** LAM ***",lam
sys.stdout.flush()

(cvode,fun_timings)=ocaml.raw_make_linalg_machine_cvode(\
    lam, field_phi0,
    "jacobian","execute_jplan",
    "rhs",
    # True, # same_nonzero_pattern for jacobian
    False, # same_nonzero_pattern for jacobian # Why doesn't "true" work here?!??
    2, # max_order
    10, # krylov_max
    )

def print_field(field,name="phi",t=None):
    for pos in range(0,11):
        probed=ocaml.probe_field(field,name,[float(pos)])
        print "T=%10s POS: %5.2f VAL: %s"%(str(t),pos,probed)

print "*** CVODE ***",cvode
sys.stdout.flush()

print_field(field_phi0,t=0.0)

for tt in range(1,100):
    target_time=tt*1.5
    ocaml.raw_cvode_advance(cvode,field_phi,target_time,-1)
    print " --- "
    print_field(field_phi,t=target_time)
