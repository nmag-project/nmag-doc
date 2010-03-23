
import logging
import nsim.logtools
import nfem.visual
import ocaml,math

from nsim import linalg_machine as nlam
import nmesh

nsim.logtools.setGlobalLogLevel("debug")

def jT_lam(name="jT",
           sigma_el=1.0, # electrical conductivity
           sigma_th=1.0, # thermal conductivity
           c_heat=1.0,   # heat capacity
           T_initial=300.0,
           mesh=None,
           ksp_tolerances={}
           ):

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

	#Thomas, not sure how you use this. Are aware that this never
	#returns 0 or 1?

        if coords[0]<-199.9 or coords[0]>199.9:
            return -2
        else:
            return -1

    def fun_T_initial(field,pos):
        print ("SAMPLING field %s AT: %s" % (field,pos))
        # return T_initial
	x,y,z=pos
	if x < 0: 
	       return 300
	else:
	       return 100
	    
        return 0.0 #should not occur
        

    
    #vector template
    mwe_V=ocaml.make_mwe("V",raw_mesh,[(0,ocaml.empty_element),(1,elem_V)],[fun_outer_region])
    #scalar template
    mwe_S=ocaml.make_mwe("S",raw_mesh,[(0,ocaml.empty_element),(1,elem_S)],[fun_outer_region])

    mwe_j=ocaml.mwe_sibling(mwe_V,"j","j/V",[("V","j")])
    mwe_phi=ocaml.mwe_sibling(mwe_S,"phi","phi/S",[("S","phi")])
    mwe_rho=ocaml.mwe_sibling(mwe_S,"rho","rho/S",[("S","rho")])
    
    mwe_T=ocaml.mwe_sibling(mwe_S,"T","T/S",[("S","T")])
    mwe_Laplace_T=ocaml.mwe_sibling(mwe_S,"Laplace_T","Laplace_T/S",[("S","Laplace_T")])
    mwe_dTdt=ocaml.mwe_sibling(mwe_S,"dTdt","dTdt/S",[("S","dTdt")])

    master_mwes_and_fields_by_name={}



    master_mwes_and_fields_by_name["T"]=(mwe_T,ocaml.raw_make_field(mwe_T,[fun_T_initial],"",""))



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
                                     "<j(k) || d/dxk phi>, k:%d" % dim),
        "op_Laplace_T":nlam.lam_operator("op_Laplace_T","Laplace_T","T",
                                         "-<d/dxk Laplace_T || d/dxk T>, k:%d" % dim,
                                         for_jacobian=True
                                         ),
        # Tricky issue: we have DNC along the contacts and NBC everywhere else.
        # XXX REPAIR: LHS field should be named "rho".
        # XXX Can our "left field equals right field" syntax already cope with that?
        "op_laplace_DNBC":nlam.lam_operator("op_laplace_DNBC","phi","phi",
                                            """ -<d/dxk phi[vol] || d/dxk phi[vol]>
                                                -<d/dxk phi[boundary=-1/1] || d/dxk phi[vol]>
                                                -<d/dxk phi[boundary=-1/1] || d/dxk phi[boundary=-1/1]>
                                                -<d/dxk phi[vol] || d/dxk phi[boundary=-1/1]>
                                                ;phi[boundary=-2/*]=phi[boundary=-2/*], k:%d""" %dim),
        "op_load_DBC":nlam.lam_operator("op_load_DBC","rho","phi",
                                        "<d/dxk rho[vol] || d/dxk phi[boundary=-2/*]>;(L||R)=(*||phi[boundary=-2/*]), k:%d" %dim),
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
                                                           ["TSTOP", "update_dTdt"],
                                                           ]),
                  "set_Laplace_T":nlam.lam_program("set_Laplace_T",
                                                   commands=[["TSTART","Laplace_T"],
                                                             ["SM*V","op_Laplace_T","v_T","v_Laplace_T"],
                                                             ["CFBOX","Laplace_T","v_Laplace_T"],
                                                             ["TSTOP","Laplace_T"]]),
                  "set_j":nlam.lam_program("set_j",
                                           commands=[["TSTART","j"],
                                                     ["SITE-WISE-IPARAMS","local_set_phi_heating",["v_phi"],[]],
                                                     ["PULL-FEM","phi","phi[boundary=-2/*]","v_phi","v_phi_boundary"],
                                                     # ["DEBUG","POS 2","v_T",1],
                                                     ["SM*V","op_load_DBC","v_phi_boundary","v_rho"],
                                                     ["SCALE","v_phi",0.0],
                                                     ["SOLVE","solve_laplace_DNBC","v_rho","v_phi"],
                                                     # ["DEBUG","POS 3","v_T",1],
                                                     ["SM*V","op_j_phi","v_phi","v_j"],
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
                                                   ["GOSUB", "update_dTdt"],
                                                   ["COLLECT","arg_dTdt","v_dTdt"],
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

mesh_filename="cons1.nmesh"
#jacobi_same_nonzero_pattern=True # Does not work yet!
jacobi_same_nonzero_pattern=False
max_order=2
krylov_max=100
max_it=1000000

my_mesh=nmesh.load(mesh_filename)
(lam,master_mwes_and_fields_by_name) = jT_lam(mesh=my_mesh)

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

heating_time=0.5
time_step=100

# Set the electrical contact potential to +/- Phi_ext at left and
# right contact:

ocaml.linalg_machine_set_iparam(lam,"Phi_ext",5.0)
ocaml.raw_cvode_advance(cvode,sundialsbuffer_final,heating_time,max_it)

# Turn off electric heating:
ocaml.linalg_machine_set_iparam(lam,"Phi_ext",0.0)
print "Known fieldnames are",master_mwes_and_fields_by_name.keys()
fields = map( lambda a: a[0],master_mwes_and_fields_by_name.values())
for field in fields:
    print "known fields are",ocaml.sys_ocamlpill_type(field)
    #Thomas, we only have one filed in master_mwes and fields -- is that right?x

for i in range(500):
    ocaml.raw_cvode_advance(cvode,sundialsbuffer_final,heating_time+i*time_step,max_it)

    (mwe,field)=master_mwes_and_fields_by_name['T']
    ocaml.linalg_machine_get_field(lam,field,"v_T")

    #Thomas, do we need to copy the data back into master_mwes_and_fields_by_name?
    nfem.visual.fields2vtkfile(master_mwes_and_fields_by_name['T'][1],'T-%04d.vtk' % i,my_mesh,format='binary')
    T0 = ocaml.probe_field(sundialsbuffer_final,"T",[0.0, 0.0, 0.0])
    print "i: %3d T: %s" % (i,repr(T0))

