"""
        (fun_name_to_mwe_and_field,
         fun_compute_H_total,
         fun_update_energies,
         make_timestepper)=ocaml.magsim_brain(["magsim-brain"], # opt petsc name
                                              self.features["do_demag"],
                                              self.mesh.raw_mesh,
                                              region_info,
                                              self.energy_factor,
                                              [], # self.intensive_params,
                                              ts_tuning_params_list,
                                              self.features["fem_only"],
                                              self.opt_lindholm_formula,
                                              self.features["ddd_use_linalg_machine"],
                                              self.extra_ccode,
                                              )
"""


def _material_template_params(mat):
    # Internal: take a material and produce the template params for py_magsim_brain.
    # This is a separate function so that we only have to change this code
    # once we get rid of that braindead idea of "simulation units" -- T.F.
    # Parameters returned:
    # ["%MAT%","%ABS_MAG%","%negJ%","%LLG_C1%","%LLG_C2%","%LLG_C3%"]

    return [mat.name,
            mat.su_Ms,
            -mat.su_exch_prefactor,
            mat.su_llg_coeff1,
            mat.su_llg_coeff2,
            mat_su_llg_normalisationfactor
            ]

# For now, we hook up the python nsim-machine via a modified
# magsim_brain interface.

def py_magsim_brain(petsc_name="nmag",
                    do_demag=True,
                    mesh=None,
                    region_info=[], # materials by region
                    energy_factor=1.0,
                    intensive_params=[],
                    ts_tuning_params=[],
                    # fem_only not provided!
                    lindholm_formula=None, # optionally, pass in Lindholm's Formula
                    # extra_ccode not provided!
                    ):

    if(mesh==None):
        raise NmagUserError, "No mesh provided!"

    if(region_info==[]):
        raise NmagUserError, "No region_info provided!"

    all_materials={}

    elem_H_demag=ocaml.make_element("H_demag",[3],3,1)
    elem_scalar=ocaml.make_element("scalar",[],3,1)

    for region in region_info:
        for mat in region:
            if(not all_materials.has_key(mat.name)):
                elem_E = ocaml.make_element("E_"+mat.name,[],3,1)
                elem_m = ocaml.make_element("m_"+mat.name,[3],3,1)
                all_materials[mat.name]=(mat,elem_E,elem_m)
    
    material_params=map(lambda xy:_material_template_params(xy[0]),all_materials.keys)

    elems_m = map(lambda region:reduce(ocaml.fuse_elements,map(lambda mat:all_materials[mat.name][2],region)),region_info)
    
    elems_E = map(lambda region:reduce(ocaml.fuse_elements,map(lambda mat:all_materials[mat.name][1],region)),region_info)
    
    elems_m_with_region_numbers=[]
    for (n,e) in enumerate(elems_m):
        elems_m_with_region_numbers.append((n,e))
        
    elems_E_with_region_numbers=[]
    for (n,e) in enumerate(elems_E):
        elems_E_with_region_numbers.append((n,e))

    elems_H_demag_with_region_numbers=[]
    for (n,e) in enumerate(elems_m):
        elems_E_with_region_numbers.append((n,elem_H_demag))

    elems_scalar_with_region_numbers=[]
    for (n,e) in enumerate(elems_m):
        elems_scalar_with_region_numbers.append((n,elem_scalar))


    mwe_m=ocaml.make_mwe("mwe_m",mesh,elems_m_with_region_numbers,None)
    mwe_E=ocaml.make_mwe("mwe_E",mesh,elems_E_with_region_numbers,None)
    mwe_H_demag=ocaml.make_mwe("mwe_H_demag",mesh,elems_H_demag_with_region_numbers,None)
    mwe_scalar=ocaml.make_mwe("mwe_scalar",mesh,elems_scalar_with_region_numbers,None)
    
    mwe_dm=ocaml.mwe_sibling(mwe_m,"dm","dm/m",map(lambda m: ("m_"+m.name,"dm_"+m.name),material_params))
    mwe_H_exch=ocaml.mwe_sibling(mwe_m,"H_exch","H_exch/m",map(lambda m: ("m_"+m.name,"H_exch_"+m.name),material_params))
    mwe_H_total=ocaml.mwe_sibling(mwe_m,"H_total","H_total/m",map(lambda m: ("m_"+m.name,"H_total_"+m.name),material_params))
    mwe_H_ext=ocaml.mwe_sibling(mwe_H_demag,"H_ext","H_ext/H_demag",[("H_demag","H_ext")])
    mwe_E_exch=ocaml.mwe_sibling(mwe_E,"E_exch","E_exch/E",map(lambda m: ("E_"+m.name,"E_exch_"+m.name),material_params))
    mwe_E_total=ocaml.mwe_sibling(mwe_E,"E_total","E_total/E",map(lambda m: ("E_"+m.name,"E_total_"+m.name),material_params))
    mwe_E_demag=ocaml.mwe_sibling(mwe_E,"E_demag","E_demag/E",map(lambda m: ("E_"+m.name,"E_demag_"+m.name),material_params))

    mwe_E_ext=ocaml.mwe_sibling(mwe_E,"E_ext","E_ext/E",map(lambda m: ("E_"+m.name,"E_ext_"+m.name),material_params))

    mwe_rho=ocaml.mwe_sibling(mwe_scalar,"rho","rho/scalar",[("scalar","rho")])
    mwe_phi=ocaml.mwe_sibling(mwe_scalar,"phi","phi/scalar",[("scalar","phi")])

    template_vars=["%MAT%","%ABS_MAG%","%negJ%","%LLG_C1%","%LLG_C2%","%LLG_C3%"]
    
    template_eq_H_total = """
H_total_%MAT%(j) <- H_ext(j) + H_demag(j) + H_exch_%MAT%(j);
"""

    template_eq_rhs = """
dm_%MAT%(i) <-
    %LLG_C1% * eps(i,j,k) * m_%MAT%(j) * H_total_%MAT%(k)
  + %LLG_C2% * eps(i,j,k) * m_%MAT%(j) * eps(k,p,q) * m_%MAT%(p) * H_total_%MAT%(q)
  + %LLG_C3% * ((-1.0) + m_%MAT%(j)*m_%MAT%(j)) * m_%MAT%(i);
"""

    template_op_H_exch = "(%negJ%)*<d/dxj H_exch_%MAT%(k)||d/dxj m_%MAT%(k)>"
    template_op_div_m = "(%ABS_MAG%)*<rho||d/dxj m_%MAT%(j)>+(%ABS_MAG%)*<rho||D/Dxj m_%MAT%(j)>"
    

    eq_H_total="%range j:3\n"+\
                ocaml.string_multifill_template_then_concat("",
                                                            template_eq_H_total,
                                                            template_vars,
                                                            material_params)
    
    eq_rhs="%range i:3, j:3, k:3, p:3, q:3;"+\
            ocaml.string_multifill_template_then_concat("",
                                                        template_eq_H_rhs,
                                                        template_vars,
                                                        material_params)
    
    str_H_exch = ocaml.string_multifill_template_then_concat(" + ",
                                                             template_op_H_exch,
                                                             template_vars,
                                                             material_params)+" ,j:3,k:3"
    
    str_div_m = ocaml.string_multifill_template_then_concat(" + ",
                                                            template_op_div_m,
                                                            template_vars,
                                                            material_params)+" ,j:3"

    # XXX TODO: provide code to update the E fields!

    master_fields_and_mwes=[]
    master_fields_and_mwes_by_name={}

    for (name,mwe) in [("m",mwe_m),("dm",mwe_dm),
                       ("rho","mwe_rho"),("phi",mwe_phi),
                       ("H_total",mwe_H_total),("H_exch",mwe_H_exch),("H_demag",mwe_H_demag),("H_ext",mwe_H_ext),
                       ("E_total",mwe_E_total),("E_exch",mwe_E_exch),("E_demag",mwe_E_demag),("E_ext",mwe_E_ext),
                       ]:
        n_m_f=(name,mwe,ocaml.raw_make_field(mwe,[],"",""))
        
        master_fields_and_mwes.append(n_m_f)
        master_fields_and_mwes_by_name[name]=n_m_f

    lam=make_linalg_machine("lam_mumag",
                            mwes=[mwe_m,mwe_dm,
                                  mwe_H_demag,mwe_H_exch,mwe_H_total,mwe_H_ext,
                                  mwe_E_demag,mwe_E_exch,mwe_E_total,mwe_E_ext,
                                  mwe_rho,mwe_phi],
                            vectors=[lam_vector(name="v_m",mwe_name="m"),
                                     lam_vector(name="v_dm",mwe_name="dm"),
                                     lam_vector(name="v_H_ext",mwe_name="H_ext"),
                                     lam_vector(name="v_H_exch",mwe_name="H_exch"),
                                     lam_vector(name="v_H_demag",mwe_name="H_demag"),
                                     lam_vector(name="v_H_total",mwe_name="H_total"),
                                     lam_vector(name="v_E_ext",mwe_name="E_ext"),
                                     lam_vector(name="v_E_exch",mwe_name="E_exch"),
                                     lam_vector(name="v_E_demag",mwe_name="E_demag"),
                                     lam_vector(name="v_E_total",mwe_name="E_total"),
                                     lam_vector(name="v_rho",mwe_name="rho"),
                                     lam_vector(name="v_rho_s",mwe_name="rho"),
                                     lam_vector(name="v_phi",mwe_name="phi"),
                                     lam_vector(name="v_phi1",mwe_name="phi"),
                                     lam_vector(name="v_phi2",mwe_name="phi"),
                                     lam_vector(name="v_phi1b",mwe_name="phi",restriction="phi[boundary=-1/*]"),
                                     lam_vector(name="v_phi2b",mwe_name="phi",restriction="phi[boundary=-1/*]"),
                                     ],
                            operators=[lam_operator("op_H_exch","H_exch","m",
                                                    str_H_exch,
                                                    for_jacobian=True),
                                       # Note: have to operate on normalized magnetization
                                   # and include magnitude factor! Also have to parametrise properly
                                   # for multiple magnetizations!
                                   lam_operator("op_div_m","rho","m",
                                                str_div_m),
                                   lam_operator("op_neg_laplace_phi","rho","phi",
                                                "<d/dxj rho || d/dxj phi>;gauge_fix:phi, j:3"
                                                # XXX add matoptions!
                                                ),
                                   lam_operator("op_grad_phi","H_demag","phi",
                                                "-<H_demag(j) || d/dxj phi>, j:3"),
                                   lam_operator("op_laplace_DBC","phi","phi",
                                                "-<d/dxj phi[vol] || d/dxj phi[vol]>;phi[boundary=-1/*]=phi[boundary=-1/*], j:3"),
                                   lam_operator("op_load_DBC","phi","phi",
                                                "<d/dxj phi[vol] || d/dxj phi[boundary=-1/*]>;(L||R)=(*||phi[boundary=-1/*]), j:3"),
                                   ],
                        bem_matrices=[lam_bem(name="BEM",
                                              is_hlib=False,
                                              #is_hlib=True,
                                              mwe_name="phi",
                                              inside_regions=[1],
                                              )],
                        ksps=[lam_ksp(name="solve_neg_laplace_phi",
                                      matrix_name="op_neg_laplace_phi"),
                              lam_ksp(name="solve_laplace_DBC",
                                      matrix_name="op_laplace_DBC"),
                              ],
                        local_operations=[lam_local("addup_H_total",
                                                    field_mwes=["H_total","H_exch","H_demag"],
                                                    equation=eq_H_total),
                                          lam_local("set_dm",
                                                    field_mwes=["m","dm","H_total"],
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
                                              commands=[["SM*V","op_H_exch","v_m","v_H_exch"],
                                                        ["CFBOX","H_exch","v_H_exch"],
                                                        ["DEBUG","Exchange","v_H_exch",3],
                                                        ]),
                                  lam_program("set_H_demag",
                                              commands=[["SM*V","op_div_m","v_m","v_rho"],
                                                        ["SCALE","v_rho",-1.0],
                                                        ["SOLVE","solve_neg_laplace_phi","v_rho","v_phi1"],
                                                        ["PULL-FEM","phi","phi[boundary=-1/*]","v_phi1","v_phi1b"],
                                                        ["DM*V","BEM","v_phi1b","v_phi2b"],
                                                        ["SM*V","op_load_DBC","v_phi2b","v_rho_s"],
                                                        ["SOLVE","solve_laplace_DBC","v_rho_s","v_phi2"],
                                                        ["PUSH-FEM","phi","phi[boundary=-1/*]","v_phi2b","v_phi2"],
                                                        ["AXPBY",1.0,"v_phi1",0.0,"v_phi"],
                                                        ["AXPBY",1.0,"v_phi2",1.0,"v_phi"],
                                                        ["SM*V","op_grad_phi","v_phi","v_H_demag"],
                                                        ["CFBOX","H_demag","v_H_demag"]
                                                        ]
                                              ),
                                  lam_program("set_H_total",
                                              commands=[
                                                         ["SITE-WISE","addup_H_total",["v_H_total","v_H_exch","v_H_demag"],[],[]]
                                                         ]),
                                  lam_program("update_H_total",
                                              commands=[["GOSUB", "set_H_demag"],
                                                        ["GOSUB", "set_H_exch"],
                                                        ["GOSUB", "set_H_total"]]),
                                  lam_program("rhs",
                                              args_fields=[["arg_dm","dm"],["arg_m","m"]], # [arg_name,mwe_name]
                                              commands=[["DISTRIB","arg_m","v_m"],
                                                        ["GOSUB", "update_H_total"],
                                                        ["SITE-WISE","set_dm",["v_m","v_dm","v_H_exch"],[],[]],
                                                        ["COLLECT","arg_dm","v_dm"]
                                                        ]),
                                  lam_program("execute_jplan",
                                              commands=[["GOSUB", "set_H_exch"],
                                                        ["JPLAN","jacobian",["v_m","v_H_exch"]]
                                                        ]),
                                  ]
                        )

    def fun_name_to_mwe_and_field(name):
        (name,mwe,field)=master_fields_and_mwes_by_name[name]
        ocaml.linalg_machine_get_field(lam,field,name)
        return (mwe,field)

    def update_energies():
        print "NOTE: update_energies() not implemented yet!"
        return None

    def compute_H_total():
        ocaml.execute_linalg_machine(lam,"update_H_total",[],[])
    
    def make_timestepper(fun_initial_M=None,fun_initial_H_ext=None):
        field_M0=None # <- Note: this is somewhat awkward! Should re-think and change that!

        if(not(fun_initial_M==None)):
            (name,mwe,field)=master_fields_and_mwes_by_name["m"]
            ocaml.set_field(field,fun_initial_M)
            ocaml.linalg_machine_set_field(lam,field,name)
            field_M0=field
            
        if(not(fun_initial_H_ext==None)):
            (name,mwe,field)=master_fields_and_mwes_by_name["H_ext"]
            ocaml.set_field(field,fun_initial_H_ext)
            ocaml.linalg_machine_set_field(lam,field,name)

        cvode=ocaml.raw_make_linalg_machine_cvode(\
            lam, field_M0,
            "jacobian","execute_jplan",
            "rhs"
            )
        
        def timestepper(iparam_vals, target_time, max_it= -1):
            # NOTE: ignoring iparam_vals here for now!
            ocaml.raw_cvode_advance(cvode,field_M0,target_time,max_it)

        return timestepper

    return (fun_name_to_mwe_and_field,
            update_energies,
            compute_H_total,
            make_timestepper)

