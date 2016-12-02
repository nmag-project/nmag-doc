
# Note: this does not include "comparisonplot":
all: run ncol plot vtk


remesh: netgen meshconvert


scalepngs:
	python ../tools/scalepng mesh.png 10


netgen:
	netgen -geofile=ellipsoid.geo -moderate -meshfiletype="Neutral Format" -meshfile=ellipsoid.neu -batchmode

meshconvert:
	nmeshimport --netgen ellipsoid.neu ellipsoid.nmesh.h5
	gzip -f ellipsoid.neu


run:
	../../../../bin/nsim ellipsoid.py --clean

ncol:
	../../../../bin/ncol ellipsoid H_ext_0 m_Py_0 > plot.dat

plot:  
	gnuplot make_plot.gnu
	convert -trim -antialias -density 150 hysteresis.eps hysteresis.png	
	python ../tools/scalepng hysteresis.png 15

vtk:
	../../../../bin/nsim ../../../../bin/nmeshpp --vtk ellipsoid.nmesh.h5 mesh.vtk
	echo "Use mayavi mesh.mv to display mesh, then grab screen"

clean:
	rm ellipsoid_dat* 


pdf:
	rst2latex --strict --stylesheet=stylesheet.tex -v doc.txt doc.tex
	pdflatex doc.tex; 

