import numpy,pyvtk
nx,ny=100,100

mfmdata = numpy.zeros((nx,ny)) #defaults to float type


def numpy2vtkgrid_file(outfilename,xmin,xmax,ymin,ymax,numpy_matrix,
                       zlevel=0.0,vtk_format='ascii',
                       vtkfilecomment="Created with numpy2vtkgrid",
                       vtkdatacomment="Scalar field"):
    assert len(numpy_matrix.shape) == 2,\
           "data needs to be 2d, but shape is '%s'" % numpy_matrix.shape

    nx,ny = numpy_matrix.shape

    xpos = numpy.linspace(xmin,xmax,nx)
    ypos = numpy.linspace(ymin,ymax,ny)

    print "point set = %d*%d=%d" % (nx,ny,nx*ny)
    
    vtk = pyvtk.VtkData(
        pyvtk.RectilinearGrid(\
            xpos.tolist(),ypos.tolist(),[zlevel]
            ),
        vtkfilecomment,
        pyvtk.PointData(
             pyvtk.Scalars( numpy_matrix.flatten().tolist(),
                            vtkdatacomment )
             ),
        format=vtk_format
        )
    vtk.tofile(outfilename,format=vtk_format)


import nmag
from nmag import SI, at

#create simulation object
sim = nmag.Simulation()

# define magnetic material
Py = nmag.MagMaterial( name="Py",
                       Ms=SI(795774,"A/m"),
                       exchange_coupling=SI(13.0e-12, "J/m")
                     )

# load mesh: the mesh dimensions are scaled by 100nm
sim.load_mesh( "../example_mfm/disk.nmesh.h5",
               [("disk", Py)],
               unit_length=SI(1e-9,"m"))


relaxed_m=nmag.get_subfield_from_h5file('disk_dat.h5','m_Py',row=0)
sim.set_m(relaxed_m)

#seems I need to do that to up-date-H_demag -> bug!
sim.advance_time(SI(1e-12,'s'))

#probe H_demag


import numpy



nx,ny=40,40
zlevel=7e-9
Hd=numpy.zeros((nx,ny))
xpos = numpy.linspace(-100e-9,100e-9,nx)
ypos = numpy.linspace(-100e-9,100e-9,ny)
for i in range(nx):
    for j in range(ny):
        #print "i=%d, j=%d, xpos=%g, ypos=%g" % (i,j,xpos[i],ypos[j])
        #nmag.ipython()
        Hd[i,j]=sim.probe_H_demag_siv([xpos[i],ypos[j],zlevel])[2]

        

numpy2vtkgrid_file('test.vtk',xpos[0]*1e9,xpos[nx-1]*1e9,ypos[0]*1e9,ypos[ny-1]*1e9,Hd)
    
     
        


 
