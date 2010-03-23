BIN=../../../../../bin
NSIM=$(BIN)/nsim 
NMESHIMPORT=$(BIN)/nmeshimport
NMESHPP=$(BIN)/nmeshpp	

sphere_in_box.nmesh.h5: sphere_in_box.msh
	$(NSIM) $(NSIMIMPORT) --gmsh sphere_in_box.msh sphere_in_box.nmesh.h5

sphere_in_box.msh: sphere_in_box.geo
	gmsh -3 -format msh1 -optimize_netgen -order 1 -o sphere_in_box.msh sphere_in_box.geo 

clean:
	rm -f sphere_in_box.msh sphere_in_box.nmesh.h5 mesh.vtk


mesh.vtk: sphere_in_box.nmesh.h5
	$(NSIM) $(NMESHPP) --vtk sphere_in_box.nmesh.h5 mesh.vtk

mayavi: mesh.vtk
	mayavi -d mesh.vtk -m SurfaceMap
