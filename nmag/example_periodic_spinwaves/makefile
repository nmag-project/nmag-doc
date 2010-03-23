all:
	make run
	make ncol
	make plot

netgen:
	netgen -geofile=periodic.geo -moderate -meshfiletype="Neutral Format" -meshfile=periodic.neu -batchmode

meshconvert:
	nmeshimport --netgen periodic.neu netgen.nmesh
	gzip -f periodic.neu

meshmirror:
	python ../../../../bin/nmeshmirror netgen.nmesh 1e-6 1e-6 -1,0,0 periodic.nmesh
run:
	../../../../bin/nsim spinwaves.py --clean
ncol:
	../../../../bin/ncol spinwaves time m_Py_0 m_Py_1 m_Py_2 > plot.txt

vtk:
	../../../../bin/nsim ../../../../bin/nmagpp --range 0 --vtk=initial spinwaves
	../../../../bin/nsim ../../../../bin/nmagpp --range -1 --vtk=last spinwaves
clean:
	rm spinwaves_dat* 


scalepngs:
	python ../tools/scalepng periodic_mesh.png 10
	python ../tools/scalepng initial_magn.png 10
	python ../tools/scalepng evolution-1.png 10
	python ../tools/scalepng evolution-2.png 10
	python ../tools/scalepng evolution-3.png 10


pdf:
	rst2latex --strict --stylesheet=stylesheet.tex -v doc.txt doc.tex
	pdflatex doc.tex;


_rst:
	echo "Only for testing: -> example_periodic_spinwaves.html"
	rst2html example_periodic_spinwaves.txt example_periodic_spinwaves.html
