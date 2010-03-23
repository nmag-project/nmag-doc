#
# (C) 2006 Dr. Thomas Fischbacher
# Exchange spring computation using nmag
#

# XXX Ad anisotropy: this is new and hot. Perhaps it would be nice to
# also provide a function to plot an energy hypersurface in 3d - so
# that we can use the very same function definition to both
# auto-generate a plot and use it in the simulation!

import nmag1 as mag, nmesh as nm, sys, math, time
import nfem # for debugging

mag.set_intensive_parameters(["H_x","H_y","H_z"], # may also set T,p,E, etc.
                             external_magnetization=["H_x","H_y","H_z"])

mag.set_features({"demag":False,"exchange":True,"timestep":True}) # XXX change usage!

mat_Dy = mag.MagMaterial("Dy",Ms=2.5,A=1.5,  # A = exchange constant
                         anisotropy=(2,lambda r: -10.0*r[2]*r[2]+6.0*r[1]*r[2]-4.0*r[2]*r[2]),
                         extra_H="H_total_Dy(0) -= mag[0]*2.5;",
                         )


mat_Fe = mag.MagMaterial("Fe",Ms=1.2,A=4.0,
                         extra_H="""H_total_Fe(0) -= mag[0]*1.2;
                         if(X(0)<0.01 || X(0)>3.99) /* pin the M field at boundaries */
                         {
                          H_total_Fe(0)=0.0;
                          H_total_Fe(1)=0.0;
                          H_total_Fe(2)=0.0;
                         }
                         """
                         )

mag.set_local_magnetic_coupling(mat_Dy,mat_Fe,-80.0)

mag.defregion("DyFe2",nm.union([nm.box([0.0],[0.5]),
                                nm.box([1.5],[2.5]),
                                nm.box([3.5],[4.0])]),
              mag_mat=[mat_Dy,mat_Fe]
              )


mag.defregion("YFe2", nm.union([nm.box([0.5],[1.5]),
                                nm.box([2.5],[3.5])]),
              mag_mat=mat_Fe # can pass either a material or a list of materials
              )

mag.set_meshing_parameters(cache_name="exchange-spring-mesh",
                           bounding_box=([-1.0],[5.0]),
                           a0=0.1,
                           max_steps=100,
                           )

mag.create_mesh()

def initial_M(dof_name,coords):
    print "DDD initial_M(dof_name=\"",dof_name,"\",coords=",coords,")\n"
    dir=dof_name[1][0]
    if dir==0:
        return 0.0
    elif dir==1:
        return math.sin(2.0*math.pi*coords[0]/8.0)
    else:
        return math.cos(2.0*math.pi*coords[0]/8.0)

mag.set_magnetization(initial_M)

def debugprint(n,stem,site,pos,val):
    print "N=",n," name=",stem," site=",site," pos=",pos," value=",val

nfem.field_entry_wise(mag.default_simulation_context.field_m,debugprint)

print mag.probe([0.1])
print mag.probe([0.2])
print mag.probe([0.3])
print mag.probe([0.4])

mag.advance_time([0.0,0.0,0.0],time=0.001)



import nfem.visual

nfem.visual.fields2vtkfile([mag.default_simulation_context.field_M],'data%05d.vtk' % 0,mesh=mag.default_simulation_context.mesh)


for i in range(1,40):
    mag.advance_time([0.0,2.5,0.0],time=0.5) # apply an external field in +y direction
    #mag.advance_time([0.0,0.0,0.0],steps=100)
    print "Total magnetization(%d): " % i, mag.integrate()
    print "Energy(%d):" %i,mag.total_energy()

    nfem.visual.fields2vtkfile([mag.default_simulation_context.field_M],'data%05d.vtk' % i,mesh=mag.default_simulation_context.mesh)
    sys.stdout.flush()

time.sleep(10000)
