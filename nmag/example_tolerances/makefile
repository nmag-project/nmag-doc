all: netgen meshconvert run ncol plots scalepngs table

remesh: netgen meshconvert

netgen:
	netgen -geofile=../example2/bar30_30_100.geo -meshfiletype="Neutral Format" -meshfile=bar30_30_100.neutral -batchmode

meshconvert:
	nmeshimport --netgen bar30_30_100.neutral bar30_30_100.nmesh.h5
	rm bar30_30_100.neutral

run:
	../../../../bin/nsim bar_tol.py --clean

ncol:
	../../../../bin/ncol bar_0.100000 time M_Py_0 M_Py_1 M_Py_2  > M_0.100000
	../../../../bin/ncol bar_0.010000 time M_Py_0 M_Py_1 M_Py_2  > M_0.010000
	../../../../bin/ncol bar_0.001000 time M_Py_0 M_Py_1 M_Py_2  > M_0.001000
	../../../../bin/ncol bar_0.000100 time M_Py_0 M_Py_1 M_Py_2  > M_0.000100	
	../../../../bin/ncol bar_0.000010 time M_Py_0 M_Py_1 M_Py_2  > M_0.000010
	../../../../bin/ncol bar_0.000001 time M_Py_0 M_Py_1 M_Py_2  > M_0.000001
	../../../../bin/ncol bar_1.000000 time M_Py_0 M_Py_1 M_Py_2  > M_1.000000

plots:
	xmgrace -nxy M_1.000000 -nxy M_0.100000 -nxy M_0.010000 -nxy M_0.001000 -nxy M_0.000100 -nxy M_0.000010 -nxy M_0.000001  -param plot1.par -hardcopy -hdevice PNG -printfile plot1.png
	xmgrace -nxy M_1.000000 -nxy M_0.100000 -nxy M_0.010000 -nxy M_0.001000 -nxy M_0.000100 -nxy M_0.000010 -nxy M_0.000001  -param plot2.par -hardcopy -hdevice PNG -printfile plot2.png
	xmgrace -nxy M_1.000000 -nxy M_0.100000 -nxy M_0.010000 -nxy M_0.001000 -nxy M_0.000100 -nxy M_0.000010 -nxy M_0.000001  -param plot3.par -hardcopy -hdevice PNG -printfile plot3.png

	mogrify -trim plot1.png
	mogrify -trim plot2.png
	mogrify -trim plot3.png

table:	
	python create_table.py

clean:
	rm -v resultsummary.txt bar_0.*_dat.ndt

scalepngs:
	python ../tools/scalepng plot1.png 15
	python ../tools/scalepng plot2.png 15
	python ../tools/scalepng plot3.png 15

