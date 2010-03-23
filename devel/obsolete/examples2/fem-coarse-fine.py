#!../bin/nmesh2

import nmesh,sys,math

def fun_T0(dof_name_indices,position):
    return math.sin(2.0*position[0]+0.7*position[1])

# For now, we use a very very simple mesh...

rings = nmesh.union([nmesh.difference(nmesh.ellipsoid([3.0,3.0],
                                                      transform=[("shift",[-2.5,0.0])]),
                                      [nmesh.ellipsoid([1.5,1.5],
                                                       transform=[("shift",[-2.5,0.0])])]),
                     nmesh.difference(nmesh.ellipsoid([3.0,3.0],
                                                      transform=[("shift",[2.5,0.0])]),
                                      [nmesh.ellipsoid([1.5,1.5],
                                                       transform=[("shift",[2.5,0.0])])])])

boxed_rings=nmesh.intersect([rings,nmesh.box([-8.0,-2.5],[8.0,2.5])])

the_mesh = nmesh.mesh(objects = [boxed_rings],
                         cache_name="rings-mesh2",
                         a0=0.5,
                         bounding_box=[[-10.0,-3.5],[10.0,3.5]],
                         neigh_force_scale = 1.,
                         density = "density=1.;",
                         initial_settling_steps = 50,
                         max_relaxation = 4,
                         max_steps=400
                         )

coarse_mesh=the_mesh.raw_mesh

(fine_mesh,coarse_fine_info)=ocaml.finer_mesh_from_coarse_mesh(3,coarse_mesh)

element_T=ocaml.make_element("T",[],2,1);

zz = ocaml.make_coarse_fine_mwe(coarse_mesh,fine_mesh,
                               "coarse_T","fine_T",
                               [(0,ocaml.empty_element),
                                (1,element_T)
                                ],
                               coarse_fine_info,
                               "")

(mwe_T_coarse,mwe_T_fine,fun_fine_coarse,fun_coarse_fine) = zz

print "OK 1\n"
sys.stdout.flush()

field_T_fine1=ocaml.raw_make_field(mwe_T_fine,
                                   [fun_T0],
                                   "" # petsc name - auto-generated if "".
                                   )

field_T_coarse2=ocaml.cofield_to_field([],fun_fine_coarse([],"",field_T_fine1))
field_T_fine3=ocaml.cofield_to_field([],fun_coarse_fine([],"",field_T_coarse2))
field_T_coarse4=ocaml.cofield_to_field([],fun_fine_coarse([],"",field_T_fine3))

ocaml.plot_scalar_field(field_T_fine1,
                        "T",
                        "/tmp/T_fine1.ps",
                        [(-1.0,[0.0,0.0,1.0]),(1.0,[0.0,1.0,0.0])],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])

ocaml.plot_scalar_field(field_T_coarse2,
                        "T",
                        "/tmp/T_coarse2.ps",
                        [(-1.0,[0.0,0.0,1.0]),(1.0,[0.0,1.0,0.0])],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])

ocaml.plot_scalar_field(field_T_fine3,
                        "T",
                        "/tmp/T_fine3.ps",
                        [(-1.0,[0.0,0.0,1.0]),(1.0,[0.0,1.0,0.0])],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])

ocaml.plot_scalar_field(field_T_coarse4,
                        "T",
                        "/tmp/T_coarse4.ps",
                        [(-1.0,[0.0,0.0,1.0]),(1.0,[0.0,1.0,0.0])],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])

print "Fine1:   ",ocaml.probe_field(field_T_fine1,"T",[0.2,0.5])
print "Coarse2: ",ocaml.probe_field(field_T_coarse2,"T",[0.2,0.5])
print "Fine3:   ",ocaml.probe_field(field_T_fine3,"T",[0.2,0.5])
print "Coarse4: ",ocaml.probe_field(field_T_coarse4,"T",[0.2,0.5])

