# This example demonstrates periodic exchange.

import os,time,sys,math
import nmesh
import ocaml
import nfem, nfem.visual

execfile("../../interface/nsim/linalg_machine.py")

#ocaml.init_hlib("/home/fangohr/build/HLib-1.3/Library/.libs/libhmatrix-1.3.so")
#ocaml.init_hlib("/home/tf/HLib-1.3/Library/.libs/libhmatrix-1.3.so")

mesh=nmesh.load_ascii("periodic/periodic.nmesh")
raw_mesh=mesh.raw_mesh

print "MESH: ",raw_mesh
sys.stdout.flush()

if ocaml.petsc_is_mpi():
    print "*** PARALLEL EXECUTION *** (not yet supported here!)"
    sys.exit(0) # DDD add parallel support!

    nr_nodes=ocaml.petsc_mpi_nr_nodes()
    nr_points=ocaml.mesh_nr_points(raw_mesh)
    z=nr_points/nr_nodes
    distrib = [int(round(z*(i+1)))-int(round(z*i)) for i in range(0,nr_nodes)]
    slack=nr_points-reduce(lambda x,y:x+y,distrib)
    distrib[0] = distrib[0] + slack
    print "*** RAW MESH %s *** DISTRIB %s ***" %(repr(raw_mesh),repr(distrib))
    ocaml.mesh_set_vertex_distribution(raw_mesh,distrib)

elem_V=ocaml.make_element("V",[3],2,1) # vector
elem_S=ocaml.make_element("S",[],2,1) # scalar element

mwe_V=ocaml.make_mwe("V",raw_mesh,[(0,ocaml.empty_element),(1,elem_V)],[]) 
mwe_S=ocaml.make_mwe("S",raw_mesh,[(0,ocaml.empty_element),(1,elem_S)],[])


mwe_M=ocaml.mwe_sibling(mwe_V,"M","M/V",[("V","M")])
mwe_dM=ocaml.mwe_sibling(mwe_V,"dM","dM/V",[("V","dM")])
mwe_H_exch=ocaml.mwe_sibling(mwe_V,"H_exch","H_exch/V",[("V","H_exch")])
mwe_H_total=ocaml.mwe_sibling(mwe_V,"H_total","H_total/V",[("V","H_total")])

mwe_rho=ocaml.mwe_sibling(mwe_S,"rho","rho/S",[("S","rho")])
mwe_phi=ocaml.mwe_sibling(mwe_S,"phi","phi/S",[("S","phi")])

def fun_M0(dof_name,dof_pos):
    ix=dof_name[1][0]
    xpos = dof_pos[0]
    ypos = dof_pos[1]
    xspos=math.sin(2*3.141592653589793*(2.0+xpos+5.0)/10.0)
    yspos=math.cos(2*3.141592653589793*(3.0+ypos+5.0)/10.0)
    r2=xspos-0.5*yspos
    alpha=math.exp(-r2*r2)
    if ix==2:
        return 0.0
    elif ix==1:
        return math.sin(alpha)
    else:
        return -math.cos(alpha)

def fun_M0_v2(dof_name,dof_pos):
    ix=dof_name[1][0]
    xpos = (dof_pos[0]+5.0)/10.0
    alpha = xpos * 2 * 3.1415926535
    if ix==0:
        return 0.0
    elif ix==1:
        return math.sin(alpha)
    else:
        return math.cos(alpha)


field_M0=ocaml.raw_make_field(mwe_M,[fun_M0],"","")
field_M=ocaml.raw_make_field(mwe_M,[],"","")
field_H_exch=ocaml.raw_make_field(mwe_H_exch,[],"","")

eq_rhs="""
%range i:3, j:3, k:3, p:3, q:3;

dM(i) <-
    1.2*eps(i,j,k)*M(j)*H_total(k)
  + (-0.03)*eps(i,j,k)*M(j)*eps(k,p,q)*M(p)*H_total(q)
  + (-1.0)*(M(j)*M(j)+(-1.0))*M(i);
""";

eq_H_total="""
%range j:3;

H_total(j) <- H_exch(j);
"""


lam=make_linalg_machine("lam_mumag",
                        mwes=[mwe_M,mwe_dM,mwe_H_exch,mwe_H_total,mwe_rho, mwe_phi],
                        vectors=[lam_vector(name="v_M",mwe_name="M"),
                                 lam_vector(name="v_dM",mwe_name="dM"),
                                 lam_vector(name="v_H_exch",mwe_name="H_exch"),
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
                                                "-2.5*<d/dxj H_exch(k)||d/dxj M(k)>; periodic:H_exch(k), j:2,k:3",
                                                #"-2.5*<d/dxj H_exch(k)||d/dxj M(k)>, j:3,k:3",
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
                                   ],
                        local_operations=[lam_local("addup_H_total",
                                                    field_mwes=["H_total","H_exch"],
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
                                                        # ["DEBUG","Exchange","v_H_exch",30],
                                                        ]),
                                  lam_program("set_H_total",
                                              commands=[
                                                         ["SITE-WISE","addup_H_total",["v_H_total","v_H_exch"],[],[]]
                                                         ]),
                                  lam_program("rhs",
                                              args_fields=[["arg_dM","dM"],["arg_M","M"],["equalized_M","M"]],
                                              commands=[
                                                        ["DISTRIB","arg_M","v_M"],
                                                        # ["DEBUG","m","v_M",30],
                                                        
                                                        ["GOSUB", "set_H_exch"],
                                                        ["GOSUB", "set_H_total"],
                                                        ["SITE-WISE","set_dM",["v_M","v_dM","v_H_exch"],[],[]],
                                                        # ["DEBUG","dm","v_dM",30],
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

# sys.exit(0) # DDD

(cvode,fun_timings)=ocaml.raw_make_linalg_machine_cvode(\
    lam, field_M0,
    "jacobian","execute_jplan",
    "rhs",
    True, # same_nonzero_pattern for jacobian
    2, # max_order
    300, # krylov_max
    )

def print_field(field,name="M",pos=[-4.0,1.0],t=None):
    probed=ocaml.probe_field(field,name,pos)
    v=probed[0][1]
    if t!=None:
        print "%8.6f  %8.6f %8.6f %8.6f" % (t,v[0],v[1],v[2])
    else:
        print "%8.6f %8.6f %8.6f" % (v[0],v[1],v[2])

print "*** CVODE ***",cvode
sys.stdout.flush()

print_field(field_M0,t=0.0)

for tt in range(1,500):
    target_time=tt*0.02
  
    ocaml.raw_cvode_advance(cvode,field_M,target_time,-1)
    print_field(field_M,t=target_time)
    ocaml.linalg_machine_get_field(lam,field_H_exch,"v_H_exch")
    ocaml.linalg_machine_get_field(lam,field_M,"v_M")
    nfem.visual.fields2vtkfile([field_M,field_H_exch],'05_h_periodic_%03d.vtk' % tt, mesh)

