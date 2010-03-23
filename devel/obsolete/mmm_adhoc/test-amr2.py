#!../bin/nmesh2

import nmesh, math,time, sys, time
import nfem2

print "TIME:",time.time()
sys.stdout.flush()

nfem2.set_default_dimension(2)
nfem2.set_default_order(1)

# Simulation parameters

sigma0=1.0
#alpha=0.01
alpha=-0.035

global m

print
"""
** Example: meshing four rings & solving the laplace equation
   for a space-dependent resistivity that depends on outer parameters

   (Presumably, we will encounter bugs at our first try)
**
"""

##### Loading the mesh #####
the_mesh = nmesh.load("ring.nmesh")

meshinfo = nmesh.tolists(the_mesh)
tmp = meshinfo[2][2]
nfem2.set_default_mesh(the_mesh)

##### Loading the mesh node positions #####
meshinfo = the_mesh.tolists()
node_position = meshinfo[0][2]
#print initial_magn



##### Making the elements... #####

# conductivity (scalar)
element_sigma     = nfem2.make_element("sigma",[]);
element_drho_by_dt= nfem2.make_element("drho_by_dt",[]);
element_phi       = nfem2.make_element("phi",[]);
element_J         = nfem2.make_element("J",[2]);
element_M         = nfem2.make_element("M",[3]);


# in make_mwe "1" is the region where I set the field
mwe_sigma      = nfem2.make_mwe("mwe_sigma",     [(1,element_sigma)])
mwe_drho_by_dt = nfem2.make_mwe("mwe_drho_by_dt",[(1,element_drho_by_dt)])
mwe_phi        = nfem2.make_mwe("mwe_phi",       [(1,element_phi)])
mwe_J          = nfem2.make_mwe("mwe_J",         [(1,element_J)])
mwe_M          = nfem2.make_mwe("mwe_M",         [(1,element_M)])


diffop_laplace=nfem2.diffop("-<d/dxj drho_by_dt|sigma|d/dxj phi>, j:2")
diffop_J=nfem2.diffop("-<J(k)|sigma|d/dxk phi>, k:2")


# Initial conductivity is spatially constant:

def fun_sigma0(dof_name_indices,position):
    return sigma0

# fill the field with the magnetization stored
# in the magnetization dictionary

def fun_fill_m(dof_name_indices,position):

#    for k, v in m.iteritems():
#        print k, v
        
    pos_str = str(position)
    x,y = pos_str[1:-2].split(",")
    ax,bx = x.split(".")
    ay,by = y.split(".")
    local_m = m[ax+"."+bx[:8]+","+ay+"."+by[:8]]
    idx = dof_name_indices[1][0]

    #print local_h[idx]
    #print len(local_m), position, local_m[idx]
    return local_m[idx]


# Dirichlet Boundary Conditions on our sample:

layer = 2.0                # thickness of the layer at the top and bottom
                           # of the sample where the potential is constant
def laplace_dbc(coords):
    if (coords[1]) > (2050.0 - layer) or (coords[1]) < (1050.0+layer) :
        # print "aa"
        return 1
    else:
        return 0

half_y_sample = 1550.0     # 1/2 length of the sample in the y direction

def laplace_dbc_values(dof_name_indices,coords):
    if(coords[1] > half_y_sample):
        return 1.0
    else:
        return -1.0    


import glob

intplFiles = glob.glob("*-intp_x.dat")
intplFiles.sort()
#print len(intplFiles)

#angles = ["88deg","80deg","70deg","60deg","50deg","47deg","85deg","75deg","65deg","55deg"]
angles = ["47deg"]
#fileIndices = [0, 11,41,76,80]
fileIndices = [76]
fields = [-1000,-900,-800,-700,-600,-500,-400,-350,-300,-250,-200,-180,-160,-140,-120,
          -100, -98, -96, -94, -92, -90, -88, -86, -84, -82, -80, -78, -76, -74, -72, -70, -68, -66, -64, -62, -60, -58, -56, -54, -52, -50, -48, -46, -44, -42, -40, -38, -36, -34, -32, -30, -28, -26, -24, -22, -20, -18, -16, -14, -12, -10,-8,-6,-4,-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100,
          120, 140, 160, 180, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000
          ]

current_time = time.time()
old_time =current_time

prematrix_laplace=nfem2.prematrix(diffop_laplace,
                                  mwe_drho_by_dt,mwe_phi,
                                  mwe_mid=mwe_sigma)

prematrix_J=nfem2.prematrix(diffop_J,
                            mwe_J,mwe_phi,
                            mwe_mid=mwe_sigma)


code_recompute_conductivity="""
double len2_J, sprod_MJ, cos2;
double len2_M=M(0)*M(0)+M(1)*M(1);

len2_J=J(0)*J(0)+J(1)*J(1);
sprod_MJ=M(0)*J(0)+M(1)*J(1);

cos2=(fabs(len2_J)<1e-8)?0.0:(sprod_MJ*sprod_MJ/(len2_M*len2_J));
sigma = sigma0 + alpha * cos2;
/* printf(\"cos2=%f sigma = %8.6f\\n \",cos2, sigma); */
"""

recompute_conductivity=nfem2.site_wise_applicator(parameter_names=["sigma0","alpha"],
                                                  # all the names of extra parameters
                                                  code=code_recompute_conductivity,
                                                  field_mwes=[mwe_sigma,mwe_J,mwe_M]
                                                  )

# computing the total current through the sample:

code_integrate_div_J="""
if(coords(1)>2048.0)
{
total_current_up+=drho_by_dt;
}
else if(coords(1)<1052.0)
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


def compute_total_current(the_field_sigma):
    laplace_solver=nfem2.laplace_solver(prematrix_laplace,
                                        dirichlet_bcs=[(-1,1,laplace_dbc)],
                                        mwe_mid=the_field_sigma)
    compute_div_J=nfem2.prematrix_applicator(prematrix_laplace,
                                             field_mid=the_field_sigma)

    cofield_drho_by_dt=nfem2.make_cofield(mwe_drho_by_dt) # make an empty field

    field_phi=laplace_solver(cofield_drho_by_dt,
                             dbc_values=laplace_dbc_values)

    cofield_div_J=compute_div_J(field_phi)
    return accumulate_J_total([0.0,0.0,0.0,0.0],cofields=[cofield_div_J])


def update_sigma(the_field_sigma,i):
    compute_J=nfem2.prematrix_applicator(prematrix_J,
                                         field_mid=the_field_sigma)
    laplace_solver=nfem2.laplace_solver(prematrix_laplace,
                                        dirichlet_bcs=[(-1,1,laplace_dbc)],
                                        mwe_mid=the_field_sigma)
    cofield_drho_by_dt=nfem2.make_cofield(mwe_drho_by_dt)
    #
    field_phi = laplace_solver(cofield_drho_by_dt,
                               dbc_values=laplace_dbc_values)
    print "Before cofield_to_field -",time.time(),"\n"
    sys.stdout.flush()
    field_J = nfem2.cofield_to_field(compute_J(field_phi))
    print "After cofield_to_field -",time.time(),"\n"

    sys.stdout.flush()


    #print "field_J contents:", nfem2.data_doftypes(field_J)

    # XXX NOTE: we should be able to make an applicator that does the
    # cofield_to_field conversion automatically!
    #
    #
    # Next, let us compute the new condictivity by site-wise operation on
    # field_J. We just overwrite field_sigma:
    #
    recompute_conductivity([sigma0,alpha],fields=[the_field_sigma,field_J,field_m])

    print "Refcount field_sigma:",ocaml.sys_refcount(the_field_sigma)
    return the_field_sigma



for angle in angles:
    intplFiles = glob.glob("*"+angle+"*-intp_x.dat")
    intplFiles.sort()

    print intplFiles.sort()

    amrDataFile = open("amr-data-"+angle+".txt","w")
    
    
#    fileList = [intplFiles[x] for x in fileIndices]
#    field = [fields[x]  for x in fileIndices]

    fileList = [intplFiles[x] for x in range(131)]
    field = [fields[x]  for x in range(131)]

    print field

    for fl, fi in zip(fileList,field):

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
#        for m in initial_magn:
#            print m

        ##### Fill magnetization dictionary #####
        m = {}

        for pos,magn in zip(node_position,initial_magn):
            pos_str = str(pos)
            x,y = pos_str[1:-2].split(",")
            ax,bx = x.split(".")
            ay,by = y.split(".")
            if "nan" in str(magn):
                #print ">>"+str(magn)+"<<<<<<<"
                m[ax+"."+bx[:8]+","+ay+"."+by[:8]] = [0.0,0.0,0.0]
            else:
#                print ax+"."+bx[:8]+","+ay+"."+by[:8], ">>"+str(magn)+"<<"
                
                m[ax+"."+bx[:8]+","+ay+"."+by[:8]] = magn



        # Later on, we will modify this field:

        field_sigma = None
        field_m = None
        ocaml.sys_check_heap()

        field_sigma=nfem2.make_field(mwe_sigma,fun_sigma0)
        
        field_m=nfem2.make_field(mwe_M, fun_fill_m)



        # Note the ignore_jumps parameter:
        # we indeed ARE interested just in the
        # surface outflow of electrical current!


        for i in range(1,5): # was: range(1,5)
            update_sigma(field_sigma,i)
            print "Forcing ocaml GC!"
            sys.stdout.flush()
            ocaml.sys_check_heap()
            #print "Iteration ",i," sigma at origin: ",nfem2.probe_field(field_sigma,[0.0,0.0])

        total_current = compute_total_current(field_sigma)
        resistance = 2./total_current[1] # resistance: potential / (downward current)
        print "Current: ",total_current
        print "Resistance: ",resistance

        amrDataFile.write(str(fi)+"  "+str(resistance)+"  "+str(total_current)+"\n")
        amrDataFile.flush()


    amrDataFile.close()

