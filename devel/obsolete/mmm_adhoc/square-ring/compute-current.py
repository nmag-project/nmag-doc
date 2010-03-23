#!../bin/nmesh2

import nmesh, math,time, sys, time
import nfem2
import nfem2.visual as nfem2_visual

print "TIME:",time.time()
sys.stdout.flush()

nfem2.set_default_dimension(3)
nfem2.set_default_order(1)

# Simulation parameters

sigma0_Py=670.
sigma0_Au=100000.
alpha=0.035

global magn_dict

print
"""
** meshing a square ring & solving the laplace equation
   for a space-dependent resistivity that depends on outer parameters
**
"""

##### Loading the mesh #####
the_mesh = nmesh.load("corrected-square-ring.nmesh")
nfem2.set_default_mesh(the_mesh)

##### Loading the mesh node positions #####
meshinfo = the_mesh.tolists()

##print "#######################################"
##print "##           ALL NODES               ##"
##print "#######################################"
##for n in range(len(meshinfo[0][2])):
##    print n, meshinfo[0][2][n]

[meshinfo_ring,
 meshinfo_contact7,
 meshinfo_contact8,
 meshinfo_contact9,
 meshinfo_contact10,
 meshinfo_contact11,
 meshinfo_contact12] = nmesh.visual.separate_parts(meshinfo, listOfParts=[(1,),(2,),(3,),(4,),(5,),(6,),(7,)])

# the node_position is used to fill the
# magnetization mwe_M and therefore
# only the nodes of the ring are used
node_position = meshinfo_ring[0][2]

##print "#######################################"
##print "##          RING NODES               ##"
##print "#######################################"
##for n in range(len(node_position)):
##    print n, node_position[n]

##### Making the elements... #####

element_sigma_Py     = nfem2.make_element("sigma_Py",[]);
element_sigma_Au     = nfem2.make_element("sigma_Au",[]);
element_drho_by_dt= nfem2.make_element("drho_by_dt",[]);
element_phi       = nfem2.make_element("phi",[]);
element_J         = nfem2.make_element("J",[3]);
element_M         = nfem2.make_element("M",[3]);


# sigma = sigma_Py for element 1 (ring) and sigma = sigma_Au for elements [2:7] (contacts)   
mwe_sigma      = nfem2.make_mwe("mwe_sigma",     [(1, element_sigma_Py)]+[(body_nr, element_sigma_Au) for body_nr in range(2,8)])
mwe_drho_by_dt = nfem2.make_mwe("mwe_drho_by_dt",[(body_nr,element_drho_by_dt) for body_nr in range(1,8)])
mwe_phi        = nfem2.make_mwe("mwe_phi",       [(body_nr,element_phi) for body_nr in range(1,8)])
mwe_J          = nfem2.make_mwe("mwe_J",         [(body_nr,element_J) for body_nr in range(1,8)])
mwe_M          = nfem2.make_mwe("mwe_M",         [(1,element_M)])


# the conductivities sum up: temporary assumption (reasonable for now)
diffop_laplace=nfem2.diffop("-<d/dxj drho_by_dt|sigma_Py|d/dxj phi>-<d/dxj drho_by_dt|sigma_Au|d/dxj phi>, j:3")
diffop_J=nfem2.diffop("-<J(k)|sigma_Py|d/dxk phi>-<J(k)|sigma_Au|d/dxk phi>, k:3")


# Initial conductivity depends on the object (ring or contacts)
def fun_sigma0(dof_name_indices,position):
##    print "fun_sigma0 -> pos:", position
    if(dof_name_indices[0]=="sigma_Au"):
        return sigma0_Au
    else:
        return sigma0_Py

# fill the field with the magnetization stored
# in the magnetization dictionary
def fun_fill_m(dof_name_indices,position):

    global magn_dict
    
    pos_str = str(position)
    x,y,z = pos_str[1:-2].split(",")
    ax,bx = x.split(".")
    ay,by = y.split(".")
    az,bz = z.split(".")

    # the string "ax+"."+bx[:8]+","+ay+"."+by[:8]+","+az+"."+bz[:8]"
    # is used to retrieve the right value from the dictionary
    
    local_m = magn_dict[ax+"."+bx[:8]+","+ay+"."+by[:8]+","+az+"."+bz[:8]]
    idx = dof_name_indices[1][0]

    return local_m[idx]

import glob
def main():

    
    intplFiles = glob.glob("*-intp_x.dat")
    intplFiles.sort()   
    #print len(intplFiles)

    angles = ["01deg"]
    #fileIndices = [0, 11,41,76,80]
    fileIndices = [10]
    fields = [-1000,-900,-800,-700,-600,-500,-400,-350,-300,-250,-200,-180,-160,-140,-120,
              -100, -98, -96, -94, -92, -90, -88, -86, -84, -82, -80, -78, -76, -74, -72, 
              -70, -68, -66, -64, -62, -60, -58, -56, -54, -52, -50, -48, -46, -44, -42,
              -40, -38, -36, -34, -32, -30, -28, -26, -24, -22, -20, -18, -16, -14, -12, 
              -10,-8,-6,-4,-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
              32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68,
              70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 
              120, 140, 160, 180, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000
              ]

    current_time = time.time()
    old_time =current_time


    for angle in angles:
        intplFiles = glob.glob("*"+angle+"*-intp_x.dat")
        intplFiles.sort()
        
        intplFiles.sort()
        
        #fileList = [intplFiles[x] for x in fileIndices]
        #field = [fields[x]  for x in fileIndices]

        print len(intplFiles)
        fileList = [intplFiles[x] for x in range(18,131)]
        field = [fields[x]  for x in range(18,131)]
        
        print field

        for fl, fi in zip(fileList,field):

            amrDataFile = open("amr-data-"+angle+"_"+str(fi)+".txt","w")

            filename = fl.split("-in")[0]

            print "**************************************", filename, "**************************************"

            current_time = time.time() - old_time
            print "TIME:", current_time
            old_time = time.time()
            
            ##### Loading the initial magnetization #####
            f = open(filename+"-intp_x.dat", "r")
            data_x= f.readlines()
            f.close()
            f = open(filename+"-intp_y.dat", "r")
            data_y= f.readlines()
            f.close()
            f = open(filename+"-intp_z.dat", "r")
            data_z= f.readlines()
            f.close()
            initial_magn = [[float(x), float(y),float(z)] for x,y,z in zip(data_x, data_y, data_z)]
            
            
            ##### Fill magnetization dictionary #####
            global magn_dict
            magn_dict = {}
    
            # loop over the nodes of the ring
            for pos,magn in zip(node_position,initial_magn):
                pos_str = str(pos)
                x,y,z = pos_str[1:-2].split(",")
                ax,bx = x.split(".")
                ay,by = y.split(".")
                az,bz = z.split(".")
        
                # the string "ax+"."+bx[:8]+","+ay+"."+by[:8]+","+az+"."+bz[:8]"
                # is used to safely save the magnetization associated to the nodes
                magn_dict[ax+"."+bx[:8]+","+ay+"."+by[:8]+","+az+"."+bz[:8]] = magn
    

        ## BISECTION ALGORITHM TO KEEP THE CURRENT CONSTANT:
            def constant_current( voltage ):




                # Dirichlet Boundary Conditions on our sample:

                def laplace_dbc(coords):

                    #                 top contact                     bottom contact
                    if (coords[1] > 990.0 and coords[0] > 420.0) or (coords[1] < -40.0) : 
                        return 1
                    else:
                        return 0

                def laplace_dbc_values(dof_name_indices,coords):
                    if(coords[1] > 500.0): # greater or smaller than 1/2height of the sample
                        return voltage
                    else:
                        return 0.0    


                prematrix_laplace=nfem2.prematrix(diffop_laplace,
                                                  mwe_drho_by_dt,mwe_phi,
                                                  mwe_mid=mwe_sigma)

                prematrix_J=nfem2.prematrix(diffop_J,
                                            mwe_J,mwe_phi,
                                            mwe_mid=mwe_sigma)


                code_recompute_conductivity="""
                if(have_sigma_Py)
                {
                 double len2_J, sprod_MJ, cos2;
                 double len2_M=M(0)*M(0)+M(1)*M(1)+M(2)*M(2);

                 len2_J=J(0)*J(0)+J(1)*J(1)+J(2)*J(2);
                 sprod_MJ=M(0)*J(0)+M(1)*J(1)+M(2)*J(2);

                 cos2=(fabs(len2_J)<1e-8)?0.0:(sprod_MJ*sprod_MJ/(len2_M*len2_J));
                 sigma_Py = sigma0_Py /( 1. + alpha * cos2);
                }
                /* printf(\"cos2=%f sigma = %8.6f\\n \",cos2, sigma); */
                """

                recompute_conductivity=nfem2.site_wise_applicator(parameter_names=["sigma0_Py","sigma0_Au","alpha"],
                                                                  # all the names of extra parameters
                                                                  code=code_recompute_conductivity,
                                                                  field_mwes=[mwe_sigma,mwe_J,mwe_M]
                                                                  )

                # computing the total current through the sample:

                code_integrate_div_J="""
                if(coords(1)>990.0 && coords(0)>420.0)
                {
                total_current_up+=drho_by_dt;
                }
                else if(coords(1)<-40.0)
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

                accumulate_J_total=nfem2.site_wise_applicator(parameter_names=["total_current_up",
                                                                               "total_current_down",
                                                                               "total_bad_current_plus",
                                                                               "total_bad_current_minus"
                                                                               ],
                                                              code=code_integrate_div_J,
                                                              position_name="coords",
                                                              cofield_mwes=[mwe_drho_by_dt])


                def update_sigma(the_field_sigma,i):
                    compute_J=nfem2.prematrix_applicator(prematrix_J,
                                                         field_mid=the_field_sigma)
                    laplace_solver=nfem2.laplace_solver(prematrix_laplace,
                                                        dirichlet_bcs=[(-1,4,laplace_dbc),(-1,6,laplace_dbc)],
                                                        mwe_mid=the_field_sigma)
                    cofield_drho_by_dt=nfem2.make_cofield(mwe_drho_by_dt)
                    #
                    field_phi = laplace_solver(cofield_drho_by_dt,
                                               dbc_values=laplace_dbc_values)

                    field_J = nfem2.cofield_to_field(compute_J(field_phi))

                    # Next, let us compute the new condictivity by site-wise operation on
                    # field_J. We just overwrite field_sigma:
                    recompute_conductivity([sigma0_Py,sigma0_Au,alpha],fields=[the_field_sigma,field_J,field_m])

                    return the_field_sigma

                def compute_total_current(the_field_sigma):

                    laplace_solver=nfem2.laplace_solver(prematrix_laplace,
                                                        # boundary conditions are defined on objects 4 and 6
                                                        dirichlet_bcs=[(-1,4,laplace_dbc),(-1,6,laplace_dbc)],
                                                        mwe_mid=the_field_sigma)
                    compute_div_J=nfem2.prematrix_applicator(prematrix_laplace,
                                                             field_mid=the_field_sigma)

                    cofield_drho_by_dt=nfem2.make_cofield(mwe_drho_by_dt)

                    field_phi=laplace_solver(cofield_drho_by_dt,
                                             dbc_values=laplace_dbc_values)

                    compute_J=nfem2.prematrix_applicator(prematrix_J,
                                                         field_mid=the_field_sigma)

                    field_J = nfem2.cofield_to_field(compute_J(field_phi))

                    cofield_div_J=compute_div_J(field_phi)

                    field_div_J = nfem2.cofield_to_field(cofield_div_J)
                    current = accumulate_J_total([0.0,0.0,0.0,0.0],cofields=[cofield_div_J])

                    return current, field_phi, the_field_sigma, field_J


                global field_m
                global field_sigma
                # Later on, we will modify this fields:
                field_sigma = None
                field_m = None

                ocaml.sys_check_heap()

                field_sigma=nfem2.make_field(mwe_sigma,fun_sigma0)
        ##        print "#######################################"
        ##        print "##          FIRST BATCH  SIGMA       ##"
        ##        print "#######################################"
        ##        nfem2.field_print_contents(field_sigma)

                field_m=nfem2.make_field(mwe_M, fun_fill_m)
        ##        print "#######################################"
        ##        print "##       FIRST BATCH  M              ##"
        ##        print "#######################################"
        ##        nfem2.field_print_contents(field_m)


                # update sigma
                for i in range(1,3): 
                    update_sigma(field_sigma,i)
                    print "Forcing ocaml GC!"
                    sys.stdout.flush()
            ##        nfem2.field_print_contents(field_sigma)
                    ocaml.sys_check_heap()


                total_current, field_phi, field_sigma, field_J = compute_total_current(field_sigma)

                nfem2_visual.fields2vtkfile( field_sigma, "sigma"+str(fi)+".vtk")
                nfem2_visual.fields2vtkfile( field_phi, "phi"+str(fi)+".vtk")
                nfem2_visual.fields2vtkfile( field_J, "current-density"+str(fi)+".vtk")

                voltage_contact7 = nfem2.probe_field(field_phi,[-20.0,630.0, 40.0])[0][1]
                voltage_contact8 = nfem2.probe_field(field_phi,[-20.0,385.0, 40.0])[0][1]

                print total_current
                resistance = voltage/total_current[1] # resistance: potential / (downward current)
                print "Current: ",total_current
                print "Resistance: ",resistance

                # ideal voltage is 150 in our units (to respect 150 microA in the experiment)
                string2file = "voltage-contacts_11-9: %f\t total current: %f\t contact7: %f\t contact8: %f\n" % (voltage, total_current[1], voltage_contact7, voltage_contact8)
                amrDataFile.write(str(fi)+"  "+str(resistance)+"  "+str(total_current)+"\n")
                amrDataFile.write(string2file)
                amrDataFile.flush()
                return total_current[1]-150.


            import scipy
            from scipy.optimize import newton
            newton (func = constant_current, x0=0.00926, tol =0.0001)
            amrDataFile.close()
   
main()
