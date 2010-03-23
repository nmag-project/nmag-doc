import os,time,sys,math
import nmesh

execfile("../interface/nsim/linalg_machine.py")

clock_time0=time.time()

objects=[nmesh.ellipsoid([3.0,3.0],[("shift",[0,0])]),
         ]

mesh = nmesh.mesh(objects=objects,
                  # a0=0.5,
                  a0=0.25,
                  bounding_box=[[-5.0,-5.0],[5.0,5.0]],
                  cache_name="linalg_machine_test.mesh",
                  )

raw_mesh=mesh.raw_mesh

print "MESH: ",raw_mesh

elem_T=ocaml.make_element("T",[],2,1)

mwe_T=ocaml.make_mwe("T",raw_mesh,[(1,elem_T)],[])
mwe_Laplace_T=ocaml.mwe_sibling(mwe_T,"Laplace_T","Laplace_T/T",[("T","Laplace_T")])
mwe_dT=ocaml.mwe_sibling(mwe_T,"dT","dT/T",[("T","dT")])

def plot_T(name,field,dof_name="T"):
    ocaml.plot_scalar_field(field,dof_name,name,
                            [(1.0,[0.0,0.0,0.0]),(3.0,[1.0,1.0,1.0])],
                            1, # plot order=1
                            1, # plot edges?
                            [-6.0,7.0,40.0,40.0])

def fun_T0(dof_name,dof_pos):
    # T=2.0+0.5*dof_pos[1] # works
    # T=2.0+math.sin(3.0*dof_pos[1]) # works
    T=2.0+0.3*dof_pos[1]*dof_pos[0] # works
    # T=2.0+math.cos(3.0*dof_pos[1]) # creates very interesting artefacts...
    print "DOF %s at %-20s: %8.4f"%(repr(dof_name),repr(dof_pos),T)
    return T

field_T0=ocaml.raw_make_field(mwe_T,[fun_T0],"","")
plot_T("/tmp/temperature_0.ps",field_T0)


field_T=ocaml.raw_make_field(mwe_T,[],"","")

lam=make_linalg_machine("lam_heat_conduction",
                        mwes=[mwe_T,mwe_dT,mwe_Laplace_T],
                        vectors=[lam_vector(name="T",mwe_name="T"),
                                 lam_vector(name="Laplace_T",mwe_name="Laplace_T"),
                                 lam_vector(name="dT",mwe_name="dT")],
                        operators=[lam_operator("op_Laplace_T","T","T","<d/dxj T||d/dxj T>, j:2")],
                        local_operations=[lam_local("set_dt",
                                                    field_mwes=["dT","Laplace_T"],
                                                    c_code="dT=-0.4*Laplace_T;"),
                                          lam_local("update_T",
                                                    field_mwes=["dT","T"],
                                                    c_code="T=T+0.001*dT;")
                                          ],
                        programs=[lam_program("compute_laplace_T",
                                              # NOTE: we could use a concept of "subroutines", where
                                              # "time_step" can use "compute_laplace_T", instead of
                                              # having to duplicate this!
                                              #
                                              # XXX NOTE TODO: a linalg_machine should also provide means to
                                              # register new scripts, making it more dynamical.
                                              # (Note: the issue of how much dynamicity one wants is tricky.
                                              # for debugguing purposes, just being able to register new
                                              # scripts on the fly should suffice...)
                                              commands=[["SM*V","op_Laplace_T","T","Laplace_T"],
                                                        ["CFBOX","Laplace_T","Laplace_T"]
                                                        ]),
                                  lam_program("time_step",
                                              commands=[["SM*V","op_Laplace_T","T","Laplace_T"],
                                                        ["CFBOX","Laplace_T","Laplace_T"],
                                                        ["SITE-WISE","set_dt",["dT","Laplace_T"],[],[]],
                                                        ["SITE-WISE","update_T",["dT","T"],[],[]],
                                                        ])]
                        )

ocaml.linalg_machine_set_field(lam,field_T0,"T")
print "*** T0 *** ", ocaml.probe_field(field_T0,"T",[0.0,0.0])


field_Laplace_T=ocaml.raw_make_field(mwe_Laplace_T,[],"","")

for n in range(1,2000):
    ocaml.linalg_machine_execute(lam,"time_step",[],[])
    if(n%50==0):
        ocaml.linalg_machine_get_field(lam,field_T,"T")
        plot_T("/tmp/temperature_%04d.ps"%n,field_T)
        print "%04d: T_total=%s"%(n,repr(ocaml.integrate_field(field_T,"T")))
        ocaml.linalg_machine_get_field(lam,field_Laplace_T,"Laplace_T")
        plot_T("/tmp/Laplace_T_%04d.ps"%n,field_Laplace_T,"Laplace_T")


ocaml.linalg_machine_get_field(lam,field_T,"T")


if False:
    for j in range(0,201):
        x=2.9*((j*0.01)-1.0)
        pos=[x,0.0]
        temp=ocaml.probe_field(field_T,"T",pos)
        print "%-20s: %s" % (pos,temp)

field_Laplace_T=ocaml.raw_make_field(mwe_Laplace_T,[],"","")
ocaml.linalg_machine_get_field(lam,field_Laplace_T,"Laplace_T")
plot_T("/tmp/Laplace_T.ps",field_Laplace_T,"Laplace_T")
