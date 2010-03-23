import logging
import nsim.logtools
import nfem.visual
import ocaml,math,sys

from nsim import linalg_machine as nlam
import nmesh

import numpy

#nsim.logtools.setGlobalLogLevel("debug")

#the constants


datafile="t_T_origin"
f=open(datafile,'w')

from nsim.si_units import SI

Ohm = SI(1,"V / A")
m = SI(1,"m")
cm = 0.01*m
nm = 1e-9*m
microOhm = 1e-6*Ohm
g = 0.001*SI(1,"kg")

rho_ = 20*microOhm*cm
sigma = 1.0/rho_
print "electric conductivity:",sigma
print "electric resistivity:",rho_

J = nsim.si_units.si.Joule
K = nsim.si_units.si.Kelvin
W = J/SI(1,"s")

C = 0.45*J/(g*K)
print "specific heat capacity",C

k = 85*W/(m*K)
print "thermal conductivity",k

#this is the density
rho = 8.7*g/cm**3
print "density",rho

print k/(rho*C)*SI(1,["K",1,"m",-2])
print 1/(rho*C*sigma)*(SI(1,["A",1,"m",-2]))**2

c_heat = C*rho #should it be the other way round?

print "Heat capacity is",c_heat

from nsim.su_units import SimulationUnits
su = SimulationUnits(scales={'A': 1e-3,'kg': 1e-27,'m': 1e-9,'s': 1e-12})

print "Simulation units are\n",su

T_initial=SI(0,"K")




#work out expected voltage:

j=1e12*SI("A/m^2")

#length
l = 1000*nm
#width
w = 50*nm
#height
h = 20*nm

A = w*h

R = rho_*l/A

print "Resistance is",R

I = j*A
print "Current is", I

U = R*I
print "Voltage is",U


dTdt = j*j/(rho*C*sigma)
print "dTdt = ",dTdt

print 70*"="


def jT_lam(name="jT",
           sigma_el=1.0, # electrical conductivity
           sigma_th=0.1, # thermal conductivity
           c_heat=1.0/10.,# heat capacity
           T_initial=300.0,
           mesh=None,
           ksp_tolerances={}
           ):

    print "Debug: c_heat = %f, sigma_el=%f, sigma_th=%f" % (c_heat,sigma_el,sigma_th)

    #sys.exit(0)
    intensive_params=["TIME","Phi_ext"]

    raw_mesh=mesh.raw_mesh
    dim = ocaml.mesh_dim(raw_mesh)

    def get_tol(name):
        if ksp_tolerances.has_key(name):
            ksp_tolerances[name]
        else:
            None

    elem_V = ocaml.make_element("V", [dim], dim, 1)
    elem_S = ocaml.make_element("S", [], dim, 1)

    def fun_outer_region(coords):

        if coords[0]<-499.9 or coords[0]>499.9:
            return -2
        else:
            return -1

    def fun_T_initial(field,pos):
        #print ("SAMPLING field %s AT: %s" % (field,pos))
        return T_initial

    def fun_T_initial_debug(field,pos):
        print ("SAMPLING field %s AT: %s" % (field,pos))
        # return T_initial
	x,y,z=pos
	if x < 0: 
	       return 300
	else:
	       return 100
	    
        return 0.0 #should not occur

    

    mwe_V=ocaml.make_mwe("V",raw_mesh,[(0,ocaml.empty_element),(1,elem_V)],[fun_outer_region])
    mwe_S=ocaml.make_mwe("S",raw_mesh,[(0,ocaml.empty_element),(1,elem_S)],[fun_outer_region])

    mwe_j=ocaml.mwe_sibling(mwe_V,"j","j/V",[("V","j")])
    mwe_phi=ocaml.mwe_sibling(mwe_S,"phi","phi/S",[("S","phi")])
    mwe_rho=ocaml.mwe_sibling(mwe_S,"rho","rho/S",[("S","rho")])
    
    mwe_T=ocaml.mwe_sibling(mwe_S,"T","T/S",[("S","T")])
    mwe_Laplace_T=ocaml.mwe_sibling(mwe_S,"Laplace_T","Laplace_T/S",[("S","Laplace_T")])
    mwe_dTdt=ocaml.mwe_sibling(mwe_S,"dTdt","dTdt/S",[("S","dTdt")])

    master_mwes_and_fields_by_name={}

    master_mwes_and_fields_by_name["T"]=(mwe_T,ocaml.raw_make_field(mwe_T,[fun_T_initial],"",""))
    master_mwes_and_fields_by_name["Laplace_T"]=(mwe_j,ocaml.raw_make_field(mwe_Laplace_T,[],"",""))
    master_mwes_and_fields_by_name["j"]=(mwe_j,ocaml.raw_make_field(mwe_j,[],"",""))
    master_mwes_and_fields_by_name["phi"]=(mwe_phi,ocaml.raw_make_field(mwe_phi,[],"",""))
    master_mwes_and_fields_by_name["rho"]=(mwe_rho,ocaml.raw_make_field(mwe_rho,[],"",""))

    # Note that at present, SITE_WISE does not work properly with restricted vectors,
    # so we will have to set the heating voltage throughout the sample, and throw away
    # the values inside the bulk then.

    ccode_heating= """
if(have_phi) {
  if(COORDS(0)<0) {phi= 0.5*Phi_ext;}
  else            {phi= -0.5*Phi_ext;}
 }
"""

    # Note ad resistive heating: integrated power generation = U*I.
    # Local power generation = E*j. As j = sigma*E, we get p = j^2/sigma,
    # so the contribution to dT/dt is j^2/sigma.

    eq_dTdt="""%%range k:%d;
dTdt <- (%.8f)*Laplace_T + (%.8f)*j(k)*j(k);
""" % (dim,
       sigma_th/c_heat,
       1.0/(c_heat*sigma_el))


    print "DDD eq_dTdt: ",eq_dTdt

    lam_mwes={"j":mwe_j, "phi":mwe_phi, "rho":mwe_rho,
              "T":mwe_T, "Laplace_T":mwe_Laplace_T,
              "dTdt":mwe_dTdt,
              }
    lam_vectors={"v_j":nlam.lam_vector(name="v_j",mwe_name="j"),
                 "v_phi":nlam.lam_vector(name="v_phi",mwe_name="phi"),
                 "v_phi_boundary":nlam.lam_vector(name="v_phi_boundary",
                                                  mwe_name="phi",
                                                  restriction="phi[boundary=-2/*]"),
                 "v_rho":nlam.lam_vector(name="v_rho",mwe_name="rho"),
                 "v_T":nlam.lam_vector(name="v_T",mwe_name="T"),
                 "v_Laplace_T":nlam.lam_vector(name="v_Laplace_T",mwe_name="Laplace_T"),
                 "v_dTdt":nlam.lam_vector(name="v_dTdt",mwe_name="dTdt"),
                 }
    
    lam_operators={\
        "op_j_phi":nlam.lam_operator("op_j_phi","j","phi",
                                     "%g*<j(k) || d/dxk phi>, k:%d" % (sigma_el,dim)),
        "op_Laplace_T":nlam.lam_operator("op_Laplace_T","Laplace_T","T",
                                         "-<d/dxk Laplace_T || d/dxk T>, k:%d" % dim
                                         ),
        # Tricky issue: we have DBC along the contacts and NBC everywhere else.
        # XXX REPAIR: LHS field should be named "rho".
        # XXX Can our "left field equals right field" syntax already cope with that?
        "op_laplace_DNBC":nlam.lam_operator("op_laplace_DNBC","phi","phi",
                                            """ -<d/dxk phi[vol] || d/dxk phi[vol]>
                                                -<d/dxk phi[boundary=-1/1] || d/dxk phi[vol]>
                                                -<d/dxk phi[boundary=-1/1] || d/dxk phi[boundary=-1/1]>
                                                -<d/dxk phi[vol] || d/dxk phi[boundary=-1/1]>
                                                ;phi[boundary=-2/*]=phi[boundary=-2/*], k:%d""" %dim),

        #"op_laplace_DNBC":nlam.lam_operator("op_laplace_DNBC","phi","phi",
        #                                    """ -<d/dxk phi || d/dxk phi>
        #                                        ;phi[boundary=-2/*]=phi[boundary=-2/*], k:%d""" %dim),


        "op_load_DBC":nlam.lam_operator("op_load_DBC","rho","phi",
                                        "<d/dxk rho || d/dxk phi[boundary=-2/*]>;(L||R)=(*||phi[boundary=-2/*]), k:%d" %dim),
        }

    lam_ksps={"solve_laplace_DNBC":nlam.lam_ksp(name="solve_laplace_DNBC",
                                                matrix_name="op_laplace_DNBC",
                                                ksp_type="gmres", pc_type="ilu",
                                                initial_guess_nonzero=True,
                                                rtol=get_tol("DNBC.rtol"),
                                                atol=get_tol("DNBC.atol"),
                                                dtol=get_tol("DNBC.dtol"),
                                                maxits=get_tol("DNBC.maxits"),
                                                )
              }

    lam_local={"local_dTdt":nlam.lam_local("local_dTdt",
                                           field_mwes=["dTdt","Laplace_T","j"],
                                           equation=eq_dTdt),
               "local_set_phi_heating":nlam.lam_local("local_set_phi_heating",
                                                      aux_args=intensive_params,
                                                      field_mwes=["phi"],
                                                      c_code=ccode_heating)
               }

    lam_jacobi={"jacobian":nlam.lam_jplan("jacobian","dTdt",
                                          ["T","Laplace_T","j"],
                                          [[],
                                           [# all contributions to d Laplace_T/dT
                                            ("operator","op_Laplace_T")],
                                           [# all contributions to dj/dT - none in this model.
                                           ]],
                                          eq_dTdt
                                          )} 

    lam_programs={"update_dTdt":nlam.lam_program("update_dTdt",
                                                 commands=[["TSTART","update_dTdt"],
                                                           ["GOSUB", "set_Laplace_T"],
                                                           ["GOSUB", "set_j"],
                                                           ["SITE-WISE-IPARAMS", "local_dTdt",["v_dTdt","v_Laplace_T","v_j"],[]],
                                                           # ["DEBUG","v_dTdt","v_dTdt",0], # seems okay!
                                                           ["TSTOP", "update_dTdt"],
                                                           ]),
                  "set_Laplace_T":nlam.lam_program("set_Laplace_T",
                                                   commands=[["TSTART","Laplace_T"],
                                                             ["SM*V","op_Laplace_T","v_T","v_Laplace_T"],
                                                             ["CFBOX","Laplace_T","v_Laplace_T"],
                                                             # ["DEBUG","v_Laplace_T","v_Laplace_T",0],
                                                             ["TSTOP","Laplace_T"]]),
                  "set_j":nlam.lam_program("set_j",
                                           commands=[["TSTART","j"],
                                                     ["SITE-WISE-IPARAMS","local_set_phi_heating",["v_phi"],[]],
                                                     ["PULL-FEM","phi","phi[boundary=-2/*]","v_phi","v_phi_boundary"],
                                                     # ["DEBUG","v_phi_boundary","v_phi_boundary",0],
                                                     ["SM*V","op_load_DBC","v_phi_boundary","v_rho"],
                                                     # ["DEBUG","v_rho","v_rho",0],
                                                     ["SCALE","v_phi",0.0],
                                                     ["SOLVE","solve_laplace_DNBC","v_rho","v_phi"],
                                                     ["PUSH-FEM","phi","phi[boundary=-2/*]","v_phi_boundary","v_phi"],
                                                     # ^ This is to add proper boundary values!
                                                     ["SM*V","op_j_phi","v_phi","v_j"],
                                                     # ["DEBUG","v_phi","v_phi",0],
                                                     # ["DEBUG","v_j","v_j",0],
                                                     ["CFBOX","j","v_j"],
                                                     ["TSTOP","j"]]),
                  "execute_jplan":nlam.lam_program("execute_jplan",
                                                   commands=[["TSTART","execute_jplan"],
                                                             ["GOSUB", "set_Laplace_T"],
                                                             ["GOSUB", "set_j"],
                                                             ["JPLAN","jacobian",["v_T","v_Laplace_T","v_j"]],
                                                             ["TSTOP","execute_jplan"]]),
                  "rhs":nlam.lam_program("rhs",
                                         args_fields=[["arg_dTdt","dTdt"],["arg_T","T"]], # [arg_name,mwe_name]
                                         commands=[["TSTART","rhs"],
                                                   ["DISTRIB","arg_T","v_T"],
                                                   # ["DEBUG","v_T","v_T",0],
                                                   ["GOSUB", "update_dTdt"],
                                                   ["COLLECT","arg_dTdt","v_dTdt"],
                                                   # ["DEBUG","v_dTdt","v_dTdt",0],
                                                   ["TSTOP","rhs"]
                                                   ]),
                  }

    lam=nlam.make_linalg_machine(
        name,
        intensive_params=intensive_params,
        mwes=lam_mwes.values(),
        vectors=lam_vectors.values(),
        operators=lam_operators.values(),
        ksps=lam_ksps.values(),
        local_operations=lam_local.values(),
        jacobi_plans=lam_jacobi.values(),
        programs=lam_programs.values()
        )

    return (lam,master_mwes_and_fields_by_name)

### Parameters (note that this all is quite ad-hoc and preliminary here):

mesh_filename="block.nmesh.h5"
#jacobi_same_nonzero_pattern=True # Does not work yet!
jacobi_same_nonzero_pattern=False
max_order=2
krylov_max=100
max_it=1000000

my_mesh=nmesh.load(mesh_filename)
(lam,master_mwes_and_fields_by_name) = jT_lam(mesh=my_mesh,
                                              sigma_el=su.of(sigma), # electrical conductivity
                                              sigma_th=su.of(k), # thermal conductivity
                                              c_heat=su.of(c_heat),# heat capacity
                                              T_initial=su.of(T_initial),
)

(mwe_T,field_T) = master_mwes_and_fields_by_name['T']

sundialsbuffer_initial = ocaml.raw_make_field(mwe_T,[],"","sundials_initial")
sundialsbuffer_final = ocaml.raw_make_field(mwe_T,[],"","sundials_final")
sundialsbuffer_starting = ocaml.raw_make_field(mwe_T,[],"","sundials_work")

ocaml.field_copy_into(field_T, sundialsbuffer_initial)
ocaml.field_copy_into(field_T, sundialsbuffer_starting)

(cvode,fun_timings)=ocaml.raw_make_linalg_machine_cvode(\
              lam, sundialsbuffer_starting,
              "jacobian","execute_jplan",
              "rhs",
              jacobi_same_nonzero_pattern,
              max_order,
              krylov_max,
              )

heatingTime = SI(2e-9,'s')
coolingTime = SI(2e-9,'s')
timeStep=SI(0.1e-9,'s')

heating_time=su.of(heatingTime)
cooling_time = su.of(coolingTime)
time_step=su.of(timeStep)

# Set the electrical contact potential to +/- Phi_ext at left and
# right contact:

print "Setting external potential"

#ocaml.linalg_machine_set_iparam(lam,"Phi_ext",5.0*100*1e12)

voltage = SI(0.2,"V")
print "Setting voltage to",voltage

ocaml.linalg_machine_set_iparam(lam,"Phi_ext",su.of(voltage))

for i in range(int(heating_time//time_step)):
    ocaml.raw_cvode_advance(cvode,sundialsbuffer_final,float(i*time_step),max_it)

    for name in master_mwes_and_fields_by_name.keys():
        (mwe,field)=master_mwes_and_fields_by_name[name]
        ocaml.linalg_machine_get_field(lam,field,"v_%s" % name)
    
    fields = map( lambda a: a[1],master_mwes_and_fields_by_name.values())
    #nfem.visual.fields2vtkfile(fields,'block-%02d.vtk' % i,my_mesh,format='binary')
    
    last_i = i
    
    j = numpy.array(ocaml.probe_field(master_mwes_and_fields_by_name['j'][1],'j',[0.,0.,0.])[0][1])

    print "j in A/m^2:",j*su.conversion_factor_of(SI("A/m^2"))

    T0 = ocaml.probe_field(sundialsbuffer_final,"T",[0.0, 0.0, 0.0])[0][1]
    print "i: %3d T: %s (su), T: %s " % (i,repr(T0),su.conversion_factor_of(SI(1,"K"))*T0)
    f.write("%g\t%g\n" % (i*time_step*su.conversion_factor_of(SI(1,'s')).value,su.conversion_factor_of(SI(1,"K")).value*T0))


print "Have written temperature at origin into t_T_origin"

print "Analytical derivative is ",dTdt

f.close()


