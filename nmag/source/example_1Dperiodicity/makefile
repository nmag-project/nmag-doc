all:
	make netgen
	make meshconvert
	make run
	make plot


netgen:
	netgen -geofile=cube.geo -moderate -meshfiletype="Neutral Format" -meshfile=cube.neu -batchmode

meshconvert:
	nmeshimport --netgen cube.neu cube1.nmesh.h5
	gzip -f cube.neu

run:
	../../../../bin/nsim periodic1.py --clean

plot:
	gnuplot make_plot.gnu
	convert -trim -antialias -density 150 periodic1_along_periodic_axis.eps periodic1_in_axis.png
	python ../tools/scalepng periodic1_in_axis.png 15
	convert -trim -antialias -density 150 periodic1_out_of_periodic_axis.eps periodic1_out_of_axis.png
	python ../tools/scalepng periodic1_out_of_axis.png 15

scalepngs:
	python ../tools/scalepng mesh.png 5

meshvis:
	nmeshpp -v ../examples/cube.nmesh.h5
	mayavi cube_nmesh.mv #take screen shot to produce cube_mesh.png

