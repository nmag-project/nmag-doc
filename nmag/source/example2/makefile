all: run1 ncol1 plot1 vtk run2 ncol2 plot2 run3 ncol3 plot3

remesh: netgen meshconvert

netgen:
	netgen -geofile=bar30_30_100.geo -meshfiletype="Neutral Format" -meshfile=bar30_30_100.neutral -batchmode

meshconvert:
	nmeshimport --netgen bar30_30_100.neutral bar30_30_100.nmesh.h5
	gzip -f bar30_30_100.neutral

run1:
	../../../../bin/nsim bar30_30_100.py --clean
run2:
	../../../../bin/nsim bar_relax.py --clean
run3:
	../../../../bin/nsim bar_relax2.py --clean

run4:
	../../../../bin/nsim read_h5.py

ncol1:
	../../../../bin/ncol bar time M_Py_0 M_Py_1 M_Py_2  > data_M.txt

plot1:
	xmgrace  -nxy data_M.txt -param xmgraceplot.par -hardcopy -hdevice PNG -printfile data_M.png
	convert -crop 789x437+3+3 data_M.png data_M.png

	xmgrace  -nxy data_M.txt -nxy oommf_M -param xmgraceplot.par -hardcopy -hdevice PNG -printfile data_M_OOMMF.png
	convert -crop 789x437+3+3 data_M_OOMMF.png data_M_OOMMF.png	

ncol2:
	ncol bar_relax time M_Py_0 M_Py_1 M_Py_2  > bar_relax_data_M.txt

plot2:
	xmgrace  -nxy bar_relax_data_M.txt -param xmgraceplot2.par -hardcopy -hdevice PNG -printfile bar_relax_data_M.png
	convert -crop 789x437+3+3 bar_relax_data_M.png bar_relax_data_M.png

ncol3:
	ncol bar_relax2 time M_Py_0 M_Py_1 M_Py_2  > bar_relax2_data_M.txt

plot3:
	xmgrace  -nxy bar_relax2_data_M.txt -param xmgraceplot2.par -hardcopy -hdevice PNG -printfile bar_relax2_data_M.png
	convert -crop 789x437+3+3 bar_relax2_data_M.png bar_relax2_data_M.png

vtk:
	../../../../bin/nmagpp --range 0 --vtk bar_initial.vtk bar
	../../../../bin/nmagpp --range 6 --vtk bar_final.vtk bar	

clean:
	rm bar_30_30_100_dat* 


scalepngs:
	python ../tools/scalepng bar30_30_100.png 15
	python ../tools/scalepng data_M.png 15
	python ../tools/scalepng data_M_OOMMF.png 15
	python ../tools/scalepng bar_initial_M.png 15
	python ../tools/scalepng bar_final_M.png 15
	python ../tools/scalepng bar_final_demag.png 15
	python ../tools/scalepng bar_relax_data_M.png 15
	python ../tools/scalepng bar_relax2_data_M.png 15

_rst:
	echo "Only for testing: -> example2.html"
	rst2html example2.txt example2.html
