all: updatefiles run


create_mesh:
	netgen -geofile=sphere1.geo -meshfiletype="Neutral Format" -meshfile=sphere1.neutral -batchmode
	nmeshimport --netgen sphere1.neutral sphere1.nmesh.h5

run:
	../../../../bin/nsim sphere1.py --clean

scalepngs:
	python ../tools/scalepng spheremesh.png 5
	python ../tools/scalepng sphere1_h5_m.png 15
	python ../tools/scalepng sphere1_h5_phi.png 15
	python ../tools/scalepng sphere1_h5_demag.png 15




