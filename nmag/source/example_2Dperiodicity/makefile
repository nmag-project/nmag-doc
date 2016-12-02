all:
	make netgen
	make meshconvert
	make run1
	make ncol1
	make run2
	make ncol2
	make plot1
	make plot2
	make vtk



netgen:
	netgen -geofile=prism.geo -fine -meshfiletype="Neutral Format" -meshfile=prism.neu -batchmode

meshconvert:
	../../../../bin/nmeshimport --netgen prism.neu prism.nmesh.h5
	gzip -f prism.neu


run1:
	../../../../bin/nsim no_periodic.py --clean

run2:
	../../../../bin/nsim periodic2.py --clean

ncol1:
	../../../../bin/ncol no_periodic time m_Py_0 m_Py_1 m_Py_2 > plot_no_periodic.dat

ncol2:
	../../../../bin/ncol periodic2 time m_Py_0 m_Py_1 m_Py_2 > plot_periodic2.dat

plot1:
	gnuplot make_plot_no_periodic.gnu
	convert -density 150 -antialias no_periodic.eps no_periodic.png
	python ../tools/scalepng no_periodic.png 15

plot2:
	gnuplot make_plot_periodic2.gnu
	convert -density 150 -antialias periodic2.eps periodic2.png
	python ../tools/scalepng periodic2.png 15

scalepngs:
	python ../tools/scalepng mesh.png 7

meshvis:
	../../../../bin/nsim ../../../../bim/nmeshpp -v prism.nmesh.h5
	mayavi cube_nmesh.mv #take screen shot to produce cube_mesh.png

vtk:
	../../../../bin/nsim ../../../../bin/nmagpp --range -1 --vtk no_periodic_relaxed.vtk no_periodic
	../../../../bin/nsim ../../../../bin/nmagpp --range -1 --vtk periodic2_relaxed.vtk periodic2

clean:
	rm periodic2_dat* no_periodic_dat*
        

_rst:
	echo "Only for testing: -> example6.html"
	rst2html example6.txt example6.html
