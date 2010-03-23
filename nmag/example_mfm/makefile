all: remesh run1 run2

remesh: netgen meshconvert


netgen:
	netgen -geofile=disk.geo -fine -meshfiletype="Neutral Format" -meshfile=disk.neu -batchmode

meshconvert:
	../../../../bin/nmeshimport --netgen disk.neu disk.nmesh.h5

vtk:
	../../../../bin/nmagpp --vtk relaxed.vtk disk
run1:
	../../../../bin/nsim disk.py --clean

run2:
	../../../../bin/nsim compute_mfm.py --clean


