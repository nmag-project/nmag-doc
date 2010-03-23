#!../bin/nmesh2

import nmesh, Numeric, math,time, sys, ocaml

# Simulation parameters

#h_x=0.8
#h_y=0.6

#h_x=0.0
#h_y=2.0

h_x=1.0
h_y=1.0

sigma0=1.0
alpha=-0.9

element_order=2

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
                                      [nmesh.ellipsoid([1.0,1.0],
                                                       transform=[("shift",[-2.5,0.0])])]),
                     nmesh.difference(nmesh.ellipsoid([3.0,3.0],
                                                      transform=[("shift",[2.5,0.0])]),
                                      [nmesh.ellipsoid([1.0,1.0],
                                                       transform=[("shift",[2.5,0.0])])])])

boxed_rings=nmesh.intersect([rings,nmesh.box([-8.0,-2.5],[8.0,2.5])])

N = 100
density = "density=1.;"

the_mesh = nmesh.mesh(objects = [boxed_rings],
                      cache_name="rings-mesh",
                      a0=0.508,
                      bounding_box=[[-10.0,-3.5],[10.0,3.5]],
                      neigh_force_scale = 1.,
                      density = density,
                      initial_settling_steps = 50,
                      max_relaxation = 4,
                      # callback=(my_function, N),
                      # max_steps=677
                      max_steps=400
                      )

print "Mesh: ",the_mesh
sys.stdout.flush()

##### Making the elements... #####

empty_element=ocaml.empty_element;

# conductivity (scalar)
element_sigma=ocaml.make_element("sigma",[],2,element_order);
# name, max indices, dim, order

# d rho/dt charge density, will serve as a Laplace equation RHS, zero
element_drho_by_dt=ocaml.make_element("drho_by_dt",[],2,element_order);

# Electrical potential
element_phi=ocaml.make_element("phi",[],2,element_order);

# Electrical current - a 2d vector
element_J=ocaml.make_element("J",[2],2,element_order);

print "element_J: ",element_J
sys.stdout.flush()


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

print "mwe_J: ",mwe_J
sys.stdout.flush()


# Our differential operators:

diffop_laplace=ocaml.make_diffop("-<d/dxj drho_by_dt|sigma|d/dxj phi>, j:2")
diffop_grad_phi=ocaml.make_diffop("<J(k)|sigma|d/dxk phi>, k:2")

# Initial conductivity is spatially constant:

def fun_sigma0(dof_name_indices,position):
    return sigma0

# Later on, we will modify this field:

field_sigma=ocaml.raw_make_field(mwe_sigma,
                                  [fun_sigma0],
                                  "" # petsc name - auto-generated if "".
                                  )

print "field_sigma: ",field_sigma
print "Sigma at origin: ",ocaml.probe_field(field_sigma,"sigma",[0.0,0.0])
sys.stdout.flush()

# Dirichlet Boundary Conditions on our sample:

def laplace_dbc(coords):
    if(abs(coords[1]) > (2.5-0.05)):
        return 1
    else:
        return 0

def laplace_dbc_value(dof_name_indices,coords):
    if(coords[1] > 0.0):
        return 0.5
    else:
        return -0.5    

# We need an empty source field for the laplace solver:

cofield_drho_by_dt=ocaml.raw_make_cofield(mwe_drho_by_dt,
                                          [], # do not sample a function
                                          "" # petsc name - auto-generated if "".
                                          )

prematrix_laplace=ocaml.raw_make_prematrix(diffop_laplace,
                                           mwe_drho_by_dt,
                                           mwe_phi,
                                           [mwe_sigma], # opt. middle field
                                           True,
                                           # ^ ignore_jumps:
                                           # don't do the calculations
                                           # for surface divergency-like matrix element contributions.
                                           # We do not need them here.
                                           )


prematrix_grad_phi=ocaml.raw_make_prematrix(diffop_grad_phi,mwe_J,mwe_phi,[mwe_sigma],True)

code_recompute_conductivity="""
double len2_J, sprod_HJ, cos2;
double len2_H=H_x*H_x+H_y*H_y;

len2_J=J(0)*J(0)+J(1)*J(1);
sprod_HJ=H_x*J(0)+H_y*J(1);

cos2=(fabs(len2_J)<1e-8)?0.0:(sprod_HJ*sprod_HJ/(len2_H*len2_J));
sigma = sigma0/(1 + alpha * cos2);
"""

last_field_phi=[None] # python hack

def update_sigma(the_field_sigma):
    compute_grad_phi=ocaml.prematrix_to_applicator([], # interface_coeffs
                                                   [the_field_sigma],
                                                   "",
                                                   prematrix_grad_phi)
    #
    #
    laplace_solver=ocaml.laplace_solver([], # if coeffs
                                        [(-1,1,laplace_dbc)], # dbcs
                                        [the_field_sigma], # mid field
                                        "", # Petsc name
                                        prematrix_laplace)
    #
    #
    result_field_phi = laplace_solver([],
                                      laplace_dbc_value,
                                      cofield_drho_by_dt);
    last_field_phi[0]=result_field_phi
    #
    #
    field_J = ocaml.cofield_to_field([],compute_grad_phi([], # no target
                                                         result_field_phi))
    #
    #
    # Next, let us compute the new condictivity by site-wise operation on
    # field_J. We just overwrite field_sigma:
    #
    recompute_conductivity=ocaml.site_wise_applicator([mwe_sigma,mwe_J], # all the fields we use
                                                      [],
                                                      ["H_x","H_y","sigma0","alpha"],
                                                      # all the names of extra parameters
                                                      "",
                                                      code_recompute_conductivity,
                                                      True
                                                      )
    #
    #
    recompute_conductivity([h_x,h_y,sigma0,alpha],[the_field_sigma,field_J],[])
    return the_field_sigma

update_sigma(field_sigma)

ocaml.plot_scalar_field(last_field_phi[0],
                        "phi",
                        "/tmp/mesh-phi-initial.ps",
                        # Simple color scheme:
                        # [(-1.2,[0.0,0.0,0.0]),(1.2,[1.0,1.0,1.0])],
                        # More elaborate color scheme - the trick is to include surface lines...
                        [(-0.5,[0.1,0.1,0.4]),
                         (-0.333333,[0.2,1.0,0.2]),
                         (-0.166666,[0.1,0.1,0.4]),
                         (0.0,[0.2,1.0,0.2]),
                         (0.166666,[0.1,0.1,0.4]),
                         (0.333333,[0.2,1.0,0.2]),
                         (0.5,[0.1,0.1,0.4]),
                         ],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])


#                        [(-0.6,[0.1,0.1,0.4]),
#                         (-0.5,[0.2,1.0,0.2]),
#                         (-0.4,[0.1,0.1,0.4]),
#                         (-0.3,[0.2,1.0,0.2]),
#                         (-0.2,[0.1,0.1,0.4]),
#                         (-0.1,[0.2,1.0,0.2]),
#                         (0.0,[0.1,0.1,0.4]),
#                         (0.1,[0.2,1.0,0.2]),
#                         (0.2,[0.1,0.1,0.4]),
#                         (0.3,[0.2,1.0,0.2]),
#                         (0.4,[0.1,0.1,0.4]),
#                         (0.5,[0.2,1.0,0.2]),
#                         (0.6,[0.1,0.1,0.4])],


for i in range(1,10):
    update_sigma(field_sigma)
    print "Iteration ",i," sigma at origin: ",ocaml.probe_field(field_sigma,"sigma",[0.0,0.0])

print "LFP: ",last_field_phi

ocaml.plot_scalar_field(last_field_phi[0],
                        "phi",
                        "/tmp/mesh-phi-final.ps",
                        [(-0.5,[0.1,0.1,0.4]),
                         (-0.333333,[0.2,1.0,0.2]),
                         (-0.166666,[0.1,0.1,0.4]),
                         (0.0,[0.2,1.0,0.2]),
                         (0.166666,[0.1,0.1,0.4]),
                         (0.333333,[0.2,1.0,0.2]),
                         (0.5,[0.1,0.1,0.4]),
                         ],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])

ocaml.plot_scalar_field(field_sigma,
                        "sigma",
                        "/tmp/mesh-sigma.ps",
                        # Simple color scheme:
                        # [(-1.2,[0.0,0.0,0.0]),(1.2,[1.0,1.0,1.0])],
                        # More elaborate color scheme - the trick is to include surface lines...
                        [(0.0,[0.3,1.0,0.3]),
                         (0.5,[1.0,1.0,0.3]),
                         (1.0,[1.0,0.3,0.3]),
                         ],
                        2, # plot order
                        1, # plot edges?
                        [-6.0,7.0,40.0,40.0])
