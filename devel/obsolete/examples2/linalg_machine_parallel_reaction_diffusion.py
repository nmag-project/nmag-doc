# A parallelizable demo script implementing a reaction/diffusion system
# ("burning sparkler")

import os,time,sys,math
import nmesh

execfile("../interface/nsim/linalg_machine.py")

objects=[nmesh.difference(nmesh.ellipsoid([3.0,3.0],[("shift",[0,0])]),
                          [nmesh.ellipsoid([2.5,2.5],[("shift",[0,0])])])
         ]

mesh = nmesh.mesh(objects=objects,
                  a0=0.25,
                  bounding_box=[[-5.0,-5.0],[5.0,5.0]],
                  cache_name="linalg_machine_rd.mesh",
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

T_max=5.0 # for plotting

elem_T=ocaml.make_element("T",[],2,1)
elem_F=ocaml.make_element("F",[],2,1)

mwe_T=ocaml.make_mwe("T",raw_mesh,[(1,elem_T)],[])
mwe_F=ocaml.make_mwe("F",raw_mesh,[(1,elem_F)],[])
mwe_Laplace_T=ocaml.mwe_sibling(mwe_T,"Laplace_T","Laplace_T/T",[("T","Laplace_T")])

def plot_scalar(name,field,dof_name,min=1.0,max=3.0):
    ocaml.plot_scalar_field(field,dof_name,name,
                            [(min,[0.2,0.2,1.0]),(max,[1.0,0.2,0.2])],
                            1, # plot order=1
                            1, # plot edges?
                            [-6.0,7.0,40.0,40.0])

def fun_T0(dof_name,dof_pos):
    alpha=math.atan2(dof_pos[0],dof_pos[1])
    T=10.0*math.exp(-alpha*alpha*4.0)
    print "T at %-20s: %8.4f"%(repr(dof_pos),T)
    return T


field_T0=ocaml.raw_make_field(mwe_T,[fun_T0],"","")
plot_scalar("/tmp/rd_T_0000.ps",field_T0,"T",min=0.0,max=T_max)

def fun_F0(dof_name,dof_pos):
    return 1.0

field_F0=ocaml.raw_make_field(mwe_T,[fun_F0],"","")
plot_scalar("/tmp/rd_F_0000.ps",field_F0,"F",min=0.0,max=1.0)


field_T=ocaml.raw_make_field(mwe_T,[],"","")
field_F=ocaml.raw_make_field(mwe_F,[],"","")

reaction_diffusion_TF="""
{
 static double
  burn_speed=0.01,
  T_activation=3.0,
  E_heat=3.0,
  diffusion=0.003;

  double dF, dT;

  dF=burn_speed*F*exp(-T_activation/T);
  dT=dF*E_heat-diffusion*Laplace_T;

  /* printf(\"T=%8.6f F=%8.6f dT=%8.6f dF=%8.6f\\n\",T,F,dT,dF); */

  F=F-dF;
  T=T+dT;
}
"""


lam=make_linalg_machine("lam_heat_conduction",
                        mwes=[mwe_T,mwe_F,mwe_Laplace_T],
                        vectors=[lam_vector(name="T",mwe_name="T"),
                                 lam_vector(name="F",mwe_name="F"),
                                 lam_vector(name="Laplace_T",mwe_name="Laplace_T"),
                                 ],
                        operators=[lam_operator("op_Laplace_T","T","T","<d/dxj T||d/dxj T>, j:2")],
                        local_operations=[
                                          lam_local("update_TF",
                                                    field_mwes=["T","F","Laplace_T"],
                                                    c_code=reaction_diffusion_TF)
                                          ],
                        programs=[lam_program("time_step",
                                              commands=[["SM*V","op_Laplace_T","T","Laplace_T"],
                                                        ["CFBOX","Laplace_T","Laplace_T"],
                                                        ["SITE-WISE","update_TF",["T","F","Laplace_T"],[],[]],
                                                        ])]
                        )

ocaml.linalg_machine_set_field(lam,field_T0,"T")
ocaml.linalg_machine_set_field(lam,field_F0,"F")

for n in range(1,5000):
    ocaml.linalg_machine_execute(lam,"time_step",[],[])
    if(n%50==0):
        ocaml.linalg_machine_get_field(lam,field_T,"T")
        plot_scalar("/tmp/rd_T_%04d.ps"%n,field_T,"T",min=0.0,max=T_max)
        ocaml.linalg_machine_get_field(lam,field_F,"F")
        plot_scalar("/tmp/rd_F_%04d.ps"%n,field_F,"F",min=0.0,max=1.0)

"""
$ time ../bin/nsim linalg_machine_parallel_reaction_diffusion.py

real    0m8.125s
user    0m6.997s
sys     0m1.057s

$ time mpirun -np 1 ../bin/nsim linalg_machine_parallel_reaction_diffusion.py

real    0m8.252s
user    0m6.976s
sys     0m1.131s

$ time mpirun -np 2 ../bin/nsim linalg_machine_parallel_reaction_diffusion.py

real    0m9.610s
user    0m7.190s
sys     0m0.952s

"""
