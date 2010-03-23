#!../bin/nmesh2

import nmesh, Numeric, math,time

h_x=1.0
h_y=1.0

sigma0=1.0
alpha=0.2

print
"""
** Example: meshing four rings & solving the laplace equation
   for a space-dependent resistivity that depends on outer parameters

   (Presumably, we will encounter bugs at our first try)
**
"""

##### Creating the mesh #####

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

N = 100
density = "density=1.;"

the_mesh = nmesh.mesh(objects = [boxed_rings],
                      cache_name="rings-mesh",
                      a0=0.3,
                      bounding_box=[[-10.0,-3.5],[10.0,3.5]],
                      neigh_force_scale = 1.,
                      density = density,
                      initial_settling_steps = 50,
                      max_relaxation = 4,
                      # callback=(my_function, N),
                      # max_steps=677
                      max_steps=200
                      )

##### Making the elements... #####

empty_element=ocaml.empty_element;

# conductivity (scalar)
element_sigma=ocaml.make_element("sigma",[],2,1);
# name, max indices, dim, order

# d rho/dt charge density, will serve as a Laplace equation RHS, zero
element_drho_by_dt=ocaml.make_element("drho_by_dt",[],2,1);

# Electrical potential
element_phi=ocaml.make_element("phi",[],2,1);

# Electrical current - a 2d vector
element_J=ocaml.make_element("J",[2],2,1);


mwe_sigma = ocaml.make_mwe("mwe_sigma",
                           the_mesh.raw_mesh,
                           [(0,empty_element),
                            (1,element_sigma)
                            ])

# Note: we may think about using mwe_sibling,
# though it's perhaps not a good idea here,
# physically speaking.

mwe_drho_by_dt = ocaml.make_mwe("mwe_drho_by_dt",
                           the_mesh.raw_mesh,
                           [(0,empty_element),
                            (1,element_drho_by_dt)
                            ])

mwe_phi = ocaml.make_mwe("mwe_phi",
                         the_mesh.raw_mesh,
                         [(0,empty_element),
                          (1,element_phi)
                          ])

mwe_J = ocaml.make_mwe("mwe_J",
                       the_mesh.raw_mesh,
                       [(0,empty_element),
                        (1,element_J)
                        ])

# Initial conductivity is spatially constant:

def fun_sigma0(dof_name_indices,position):
    return sigma0

field_sigma0=ocaml.raw_make_field(mwe_sigma,
                                  [fun_sigma0],
                                  "" # petsc name - auto-generated if "".
                                  )

ocaml.plot_scalar_field(field_sigma0,
                        "sigma",
                        "/tmp/mesh-sigma0.ps",
                        [(0.0,[0.0,0.0,0.0]),(2.0,[1.0,1.0,1.0])],
                        1, # plot order=1
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])


cofield_drho_by_dt=ocaml.raw_make_cofield(mwe_drho_by_dt,
                                          [],
                                          "" # petsc name - auto-generated if "".
                                          )

diffop_laplace=ocaml.make_diffop("-<d/dxj drho_by_dt|sigma|d/dxj phi>, j:2")


diffop_grad_phi=ocaml.make_diffop("<J(k)|sigma|d/dxk phi>, k:2")

prematrix_laplace=ocaml.raw_make_prematrix(diffop_laplace,
                                           mwe_drho_by_dt,
                                           mwe_phi,
                                           [mwe_sigma],
                                           1
                                           # ignore_jumps:
                                           # don't do the calculations
                                           # for surface divergency-like matrix element contributions.
                                           # We do not need them here.
                                           )

prematrix_grad_phi=ocaml.raw_make_prematrix(diffop_grad_phi,
                                            mwe_J,
                                            mwe_phi,
                                            [mwe_sigma],
                                            True
                                            )

#print "Prematrix pill type: ",ocaml.sys_ocamlpill_type(prematrix_grad_phi),"\n"


compute_grad_phi=ocaml.prematrix_to_applicator([], # interface_coeffs
                                               [field_sigma0],
                                               "grad_phi", # petsc name
                                               prematrix_grad_phi)


def laplace_dbc(coords):
    if(abs(coords[1]) > (2.5-0.05)):
        return 1
    else:
        return 0


def laplace_dbc_value(dof,coords):
    print "Computing DBC value for pos=",coords
    if(coords[1] > 0.0):
        return 1.0
    else:
        return -1.0    

laplace_solver=ocaml.laplace_solver([],
                                    # interface coefficients for boundary
                                    # jumps in the generalized laplace
                                    # operator
                                    [(-1,1,laplace_dbc)],
                                    # Identify Dirichlet Boundary Conditions:
                                    # That part of the boundary between
                                    # region 0 and region 1 which has
                                    # y-coord +/- 2.5
                                    [field_sigma0],
                                    # We do have an "operator middle field"
                                    # (otherwise, we would pass []),
                                    # and it is field_sigma0
                                    "Phi_E_Solver", # Petsc name
                                    prematrix_laplace)


result_field_phi = laplace_solver([],
                                  # we could provide [target_field_phi] here
                                  # to have the system write into that vector,
                                  # or omit it (like here) to have the system
                                  # allocate a new target vector.
                                  laplace_dbc_value,
                                  cofield_drho_by_dt);

print "Have result_field_phi:", result_field_phi

print "Phi at origin: ",ocaml.probe_field(result_field_phi,"phi",[0.0,0.0])

for n in range(0,11):
    ypos = -0.5+n*0.1
    print "Phi at ",ypos,": ",ocaml.probe_field(result_field_phi,"phi",[0.0,ypos])

print "Phi_upper: ",ocaml.probe_field(result_field_phi,"phi",[3.0,2.45])
print "Phi_lower: ",ocaml.probe_field(result_field_phi,"phi",[-3.0,2.45])

ocaml.plot_scalar_field(result_field_phi,
                        "phi",
                        "/tmp/mesh-phi.ps",
                        # Simple color scheme:
                        # [(-1.2,[0.0,0.0,0.0]),(1.2,[1.0,1.0,1.0])],
                        # More elaborate color scheme - the trick is to include surface lines...
                        [(-1.2,[0.4,0.8,0.4]),
                         (-1.0,[0.4,0.4,0.8]),
                         (-0.8,[0.4,0.85,0.4]),
                         (-0.6,[0.4,0.4,0.85]),
                         (-0.4,[0.4,0.9,0.4]),
                         (-0.2,[0.4,0.4,0.9]),
                         (0.0,[0.4,0.95,0.4]),
                         (0.2,[0.4,0.4,0.95]),
                         (0.4,[0.4,1.0,0.4]),
                         (0.6,[0.4,0.4,1.0]),
                         (0.8,[0.6,1.0,0.6]),
                         (1.0,[0.6,0.6,1.0]),
                         (1.2,[0.8,1.0,0.8])],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])


# Next, we have to compute the electrical current from this (grad phi).
# Again, we could tell the system to use a target field to write data
# to, but we choose instead to just let it allocate a new one...

field_J = ocaml.cofield_to_field([],compute_grad_phi([], # no target
                                                     result_field_phi))

# Next, let us compute the new condictivity by site-wise operation on
# field_J. We just overwrite field_sigma0:


code_recompute_conductivity="""
double len2_J, sprod_HJ, cos2;
double len2_H=H_x*H_x+H_y*H_y;

len2_J=J(0)*J(0)+J(1)*J(1);
sprod_HJ=H_x*J(0)+H_y*J(1);

cos2=(fabs(len2_J)<1e-8)?0.0:(sprod_HJ*sprod_HJ/(len2_H*len2_J));
printf(\"cos2=%f\\n\",cos2);
sigma = sigma0 - alpha * cos2; /* In this model computation, conductivity goes DOWN due to magnetic effects */
printf(\"Sigma: %f l2j: %f l2h: %f\\n\",sigma, len2_J, len2_H);fflush(stdout);

"""

recompute_conductivity=ocaml.site_wise_applicator([mwe_sigma,mwe_J], # all the fields we use
                                                  [], # cofields
                                                  ["H_x","H_y","sigma0","alpha"], # all the names of extra parameters
                                                  "",
                                                  code_recompute_conductivity,
                                                  True
                                                  )

recompute_conductivity([h_x,h_y,sigma0,alpha],[field_sigma0,field_J],[])

print "*** sigma at origin: ",ocaml.probe_field(field_sigma0,"sigma",[0.0,0.0])

ocaml.plot_scalar_field(field_sigma0,
                        "sigma",
                        "/tmp/mesh-sigma-1.ps",
                        # Simple color scheme:
                        # [(-1.2,[0.0,0.0,0.0]),(1.2,[1.0,1.0,1.0])],
                        # More elaborate color scheme - the trick is to include surface lines...
                        # [(0.0,[0.4,0.4,0.4]),
                         #(1.0,[0.4,0.8,0.4]),
                         #(1.001,[0.4,0.4,0.8]),
                         #(1.002,[0.4,0.85,0.4]),
                         #(1.003,[0.4,0.4,0.85]),
                         #(1.004,[0.4,0.9,0.4]),
                         #(1.005,[0.4,0.4,0.9]),
                         #(1.006,[0.4,0.95,0.4]),
                         #(1.007,[0.4,0.4,0.95]),
                         #(1.008,[0.4,1.0,0.4]),
                         #(1.009,[0.4,0.4,1.0]),
                         #(1.010,[0.6,1.0,0.6]),
                         #(1.011,[0.6,0.6,1.0]),
                         # (2.0,[0.8,0.8,0.8])
                         # ],
                        [(0.0,[0.4,0.4,0.4]),
                         (1.0,[0.4,1.0,0.4]),
                         (1.005,[1.0,1.0,0.4]),
                         (1.010,[1.0,0.4,0.4]),
                         (2.0,[0.8,0.8,0.8])
                         ],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])



# We need a new laplace solver that takes into account our new conductivity...

laplace_solver=ocaml.laplace_solver([],
                                    # interface coefficients for boundary
                                    # jumps in the generalized laplace
                                    # operator
                                    [(-1,1,laplace_dbc)],
                                    # Identify Dirichlet Boundary Conditions:
                                    # That part of the boundary between
                                    # region 0 and region 1 which has
                                    # y-coord 1500
                                    [field_sigma0],
                                    # We do have an "operator middle field"
                                    # (otherwise, we would pass []),
                                    # and it is field_sigma0
                                    "Phi_E_Solver_sigma1", # Petsc name
                                    prematrix_laplace)

result_field_phi1 = laplace_solver([],
                                   laplace_dbc_value,
                                   cofield_drho_by_dt);

field_J1 = compute_grad_phi([], # no target
                            result_field_phi1)

ocaml.plot_scalar_field(result_field_phi1,
                        "phi",
                        "/tmp/mesh-phi-1.ps",
                        # Simple color scheme:
                        # [(-1.2,[0.0,0.0,0.0]),(1.2,[1.0,1.0,1.0])],
                        # More elaborate color scheme - the trick is to include surface lines...
                        [(-1.2,[0.4,0.8,0.4]),
                         (-1.0,[0.4,0.4,0.8]),
                         (-0.8,[0.4,0.85,0.4]),
                         (-0.6,[0.4,0.4,0.85]),
                         (-0.4,[0.4,0.9,0.4]),
                         (-0.2,[0.4,0.4,0.9]),
                         (0.0,[0.4,0.95,0.4]),
                         (0.2,[0.4,0.4,0.95]),
                         (0.4,[0.4,1.0,0.4]),
                         (0.6,[0.4,0.4,1.0]),
                         (0.8,[0.6,1.0,0.6]),
                         (1.0,[0.6,0.6,1.0]),
                         (1.2,[0.8,1.0,0.8])],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])



# We now have the final current distribution.
# All that is left to do is to integrate the current density
# over one of the bonding contacts to get total current, and hence
# total resistivity.

# time.sleep(100000)
