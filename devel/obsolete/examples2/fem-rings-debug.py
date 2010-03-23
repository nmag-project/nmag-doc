#!../bin/nmesh2

import nmesh, Numeric, math,time, sys
import nfem

nfem.set_default_dimension(2)
nfem.set_default_order(1)
#nfem.set_default_order(2)

# Simulation parameters

h_x=1.6
h_y=1.2

sigma0=1.0
alpha=0.01


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

nfem.set_default_mesh(the_mesh)

##### Making the elements... #####

empty_element=ocaml.empty_element;

# conductivity (scalar)
element_sigma     = nfem.make_element("sigma",[]);
element_drho_by_dt= nfem.make_element("drho_by_dt",[]);
element_phi       = nfem.make_element("phi",[]);
element_J         = nfem.make_element("J",[2]);

mwe_sigma      = nfem.make_mwe("mwe_sigma",     [(1,element_sigma)])
mwe_drho_by_dt = nfem.make_mwe("mwe_drho_by_dt",[(1,element_drho_by_dt)])
mwe_phi        = nfem.make_mwe("mwe_phi",       [(1,element_phi)])
mwe_J          = nfem.make_mwe("mwe_J",         [(1,element_J)])

diffop_laplace=nfem.diffop("-<d/dxj drho_by_dt|sigma|d/dxj phi>, j:2")
diffop_J=nfem.diffop("<J(k)|sigma|d/dxk phi>, k:2")
diffop_div_J=nfem.diffop("<drho_by_dt||d/dxk J(k)>, k:2")
# ^ We need this last one to compute total current by summing over
# all the surface jumps in outflowing current.

# Initial conductivity is spatially constant:

def fun_sigma0(dof_name_indices,position):
    # For debugging, we use a somewhat more crazy conductivity which goes down
    # at the outer regions...
    r2=position[0]*position[0]+position[1]*position[1]
    factor=1.0/(1+r2)
    return (sigma0*factor)

# Later on, we will modify this field:

field_sigma=nfem.make_field(mwe_sigma,fun_sigma0)

# Dirichlet Boundary Conditions on our sample:

def laplace_dbc(coords):
    if(abs(coords[1]) > (2.5-0.05)):
        return 1
    else:
        return 0

def laplace_dbc_values(dof_name_indices,coords):
    if(coords[1] > 0.0):
        return 1.0
    else:
        return -1.0    


cofield_drho_by_dt=nfem.make_cofield(mwe_drho_by_dt)

prematrix_laplace=nfem.prematrix(diffop_laplace,
                                  mwe_drho_by_dt,mwe_phi,
                                  mwe_mid=mwe_sigma)

prematrix_J=nfem.prematrix(diffop_J,
                                   mwe_J,mwe_phi,
                                   mwe_mid=mwe_sigma)


prematrix_div_J=nfem.prematrix(diffop_div_J,
                                mwe_drho_by_dt,mwe_J,
                                ignore_jumps=False)
# Note the ignore_jumps parameter:
# we indeed ARE interested just in the
# surface outflow of electrical current!


code_recompute_conductivity="""
double len2_J, sprod_HJ, cos2;
double len2_H=H_x*H_x+H_y*H_y;

len2_J=J(0)*J(0)+J(1)*J(1);
sprod_HJ=H_x*J(0)+H_y*J(1);

cos2=(fabs(len2_J)<1e-8)?0.0:(sprod_HJ*sprod_HJ/(len2_H*len2_J));
sigma = sigma0 + alpha * cos2;
"""


compute_J=nfem.prematrix_applicator(prematrix_J,
                                     mwe_mid=field_sigma)

compute_laplace=nfem.prematrix_applicator(prematrix_laplace,mwe_mid=field_sigma)

laplace_solver=nfem.laplace_solver(prematrix_laplace,
                                    dirichlet_bcs=[(-1,1,laplace_dbc)],
                                    mwe_mid=field_sigma)

field_phi = laplace_solver(cofield_drho_by_dt,
                           dbc_values=laplace_dbc_values)

nfem.plot_scalar_field(field_phi,
                        "phi","/tmp/plot-phi.ps",
                        plot_edges=False,
                        color_scheme=[(-1.2,[0.1,0.1,0.1]),
                                      (-1.1,[0.9,0.9,0.9]),
                                      (-1.0,[0.1,0.1,0.1]),
                                      (-0.9,[0.9,0.9,0.9]),
                                      (-0.8,[0.1,0.1,0.1]),
                                      (-0.7,[0.9,0.9,0.9]),
                                      (-0.6,[0.1,0.1,0.1]),
                                      (-0.5,[0.9,0.9,0.9]),
                                      (-0.4,[0.1,0.1,0.1]),
                                      (-0.3,[0.9,0.9,0.9]),
                                      (-0.2,[0.1,0.1,0.1]),
                                      (-0.1,[0.9,0.9,0.9]),
                                      (0.0,[0.1,0.1,0.1]),
                                      (0.1,[0.9,0.9,0.9]),
                                      (0.2,[0.1,0.1,0.1]),
                                      (0.3,[0.9,0.9,0.9]),
                                      (0.4,[0.1,0.1,0.1]),
                                      (0.5,[0.9,0.9,0.9]),
                                      (0.6,[0.1,0.1,0.1]),
                                      (0.7,[0.9,0.9,0.9]),
                                      (0.8,[0.1,0.1,0.1]),
                                      (0.9,[0.9,0.9,0.9]),
                                      (1.0,[0.1,0.1,0.1]),
                                      (1.1,[0.9,0.9,0.9])])

cofield_div_J =compute_laplace(field_phi)

# This cofield should be zero everywhere except at the boundaries. Check that!

code_check_div_J="""
printf(\"%s div J at %8.4f %8.4f = %12.8f\\n\",fabs(coords(1))>2.45?\"***\":\"\",coords(0),coords(1),drho_by_dt);
fflush(stdout);
"""

code_integrate_div_J="""
printf(\"div J contrib at (%6.4f %6.4f): %f\\n\",coords(0),coords(1),drho_by_dt);fflush(stdout); /* DDD */

if(coords(1)>2.45)
{
  printf(\"Outflow Y-TOP contrib at (%6.4f %6.4f): %f\\n\",coords(0),coords(1),drho_by_dt);fflush(stdout); /* DDD */
  total_current_up+=drho_by_dt;
}
else if(coords(1)<-2.45)
{
  total_current_down+=drho_by_dt;
}
else
{
  double j=drho_by_dt;

  if(j>0)total_bad_current_plus+=j;
  else   total_bad_current_minus+=j;
}
"""

#check_div_J=nfem.site_wise_applicator(parameter_names=[],
#                                       # all the names of extra parameters
#                                       code=code_check_div_J,
#                                       cofields=[cofield_div_J],
#                                       position_name="coords"
#                                       )
#
#check_div_J([])

integrate_div_J=nfem.site_wise_applicator(parameter_names=["total_current_up",
                                                            "total_current_down",
                                                            "total_bad_current_plus",
                                                            "total_bad_current_minus"
                                                            ],
                                           code=code_integrate_div_J,
                                           position_name="coords",
                                           cofields=[cofield_div_J])



print "div J integrated: ",integrate_div_J([0.0,0.0,0.0,0.0])

sys.exit()

# This example shows that we indeed are on the right track and the
# artefacts I encuntered yesterday were just a product of working too
# long... The important thing to remember is:
# DO NOT COMPUTE div(sigma*grad(phi)), compute Laplace(phi)
# (with modified laplace operator) to get div j. The div-grad-thingy will
# do strange things when taking the gradient of phi at the boundary.


###### END ######



field_J = nfem.cofield_to_field(compute_J(field_phi))

recompute_conductivity=nfem.site_wise_applicator(parameter_names=["H_x","H_y","sigma0","alpha"],
                                                  # all the names of extra parameters
                                                  code=code_recompute_conductivity,
                                                  fields=[field_sigma,field_J]
                                                  )
#
#
recompute_conductivity([h_x,h_y,sigma0,alpha])

#for i in range(1,2): # was: range(1,5)
#    update_sigma(field_sigma)
#    print "Iteration ",i," sigma at origin: ",nfem.probe_field(field_sigma,"sigma",[0.0,0.0])

# === XXX The part below has to be re-done. For now, we just test whether it works! ===

# computing the total current through the sample:

code_integrate_div_J="""
printf(\"div J contrib at (%6.4f %6.4f): %f\\n\",coords(0),coords(1),drho_by_dt);fflush(stdout); /* DDD */

if(coords(1)>2.45)
{
  printf(\"Outflow Y-TOP contrib at (%6.4f %6.4f): %f\\n\",coords(0),coords(1),drho_by_dt);fflush(stdout); /* DDD */
  total_current_up+=drho_by_dt;
}
else if(coords(1)<-2.45)
{
  total_current_down+=drho_by_dt;
}
else
{
  double j=drho_by_dt;

  if(j>0)total_bad_current_plus+=j;
  else   total_bad_current_minus+=j;
}
"""

compute_div_J=nfem.prematrix_applicator(prematrix_div_J,
                                         interface_coeffs=[(-1,1,1.0)])

# These definitions are to be regarded as internal to compute_total_current:
_cofield_div_J=nfem.make_cofield(mwe_drho_by_dt)

_accumulate_J_total=nfem.site_wise_applicator(parameter_names=["total_current_up",
                                                                "total_current_down",
                                                                "total_bad_current_plus",
                                                                "total_bad_current_minus"
                                                                ],
                                               code=code_integrate_div_J,
                                               position_name="coords",
                                               cofields=[_cofield_div_J])

def compute_total_current(the_field_J):
    compute_div_J(the_field_J,target=_cofield_div_J)
    print "Computed div J!"
    sys.stdout.flush()
    # ^ Works up to here!
    return _accumulate_J_total([0.0,0.0,0.0,0.0])

compute_J=nfem.prematrix_applicator(prematrix_J,
                                     mwe_mid=field_sigma)

laplace_solver=nfem.laplace_solver(prematrix_laplace,
                                    dirichlet_bcs=[(-1,1,laplace_dbc)],
                                    mwe_mid=field_sigma)

field_phi = laplace_solver(cofield_drho_by_dt,
                           dbc_values=laplace_dbc_values)

field_J = nfem.cofield_to_field(compute_J(field_phi))

print "Total current: ",compute_total_current(field_J)

nfem.plot_scalar_field(nfem.cofield_to_field(_cofield_div_J),
                        "drho_by_dt","/tmp/plot-div-j.ps",
                        color_scheme=[(-2.0,[0.2,0.2,1.0]),
                                      (-0.2,[0.2,1.0,1.0]),
                                      (0.0,[0.4,0.4,0.4]),
                                      (0.2,[1.0,1.0,0.0]),
                                      (2.0,[1.0,0.2,0.2])]);


# 
# Results from a timed test run with cached mesh
# 
# element order=1
# ===============
# 
# Iteration  1  sigma at origin:  1.00359354396
# Iteration  2  sigma at origin:  1.00361510401
# Iteration  3  sigma at origin:  1.00361511452
# Iteration  4  sigma at origin:  1.00361511457
# 
# real    0m13.726s
# user    0m11.077s
# sys     0m2.612s
# 
# element order=2
# ===============
# 
# Iteration  1  sigma at origin:  1.00360032758
# Iteration  2  sigma at origin:  1.00361542514
# Iteration  3  sigma at origin:  1.00361543821
# Iteration  4  sigma at origin:  1.00361543833
# 
# real    4m15.084s
# user    3m10.796s
# sys     1m4.053s
#
