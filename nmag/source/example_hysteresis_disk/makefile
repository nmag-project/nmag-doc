all: run ncol plot comparisonplot vtk

remesh: netgen meshconvert


netgen:
	netgen -geofile=the_nanodot.geo -fine -meshfiletype="Neutral Format" -meshfile=the_nanodot.neu -batchmode

meshconvert:
	../../../../bin/nmeshimport --netgen the_nanodot.neu nanodot1.nmesh.h5
	gzip -f the_nanodot.neu

run:
	../../../../bin/nsim nanodot1.py --clean

ncol:
	../../../../bin/ncol nanodot1 H_ext_0 m_Py_0 > plot.dat
	head -n 5 plot.dat > plot_head.txt

plot:
	gnuplot make_plot.gnu
	convert -trim -antialias -density 150 nanodot_hyst.eps nanodot_hyst.png 	
	python ../tools/scalepng nanodot_hyst.png 15

comparisonplot:
	gnuplot make_comparison_plot.gnu
	convert -trim -antialias -density 150 nanodot_comparison_hyst.eps nanodot_comparison_hyst.png 
	python ../tools/scalepng nanodot_comparison_hyst.png 15


meshvis:
	../../../../bin/nsim ../../../../bin/nmeshpp -Q ../example_hysteresis/nanodot1.nmesh.h5

vtk:
	../../../../bin/nsim ../../../../bin/nmagpp --range 0 --vtk nanodot_initial.vtk nanodot1
	../../../../bin/nsim ../../../../bin/nmagpp --range -1 --vtk nanodot_final.vtk nanodot1

clean:
	rm nanodot_dat* 

scalepngs:
	python ../tools/scalepng nanodot1_mesh.png 15


_rst:
	echo "Only for testing: -> example3.html"
	rst2html example3.txt example3.html

pdf:
	rst2latex --strict --stylesheet=stylesheet.tex -v doc.txt doc.tex
	pdflatex doc.tex;

