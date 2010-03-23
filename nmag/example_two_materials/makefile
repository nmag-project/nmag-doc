netgen:
	netgen -geofile=two_cubes.geo -meshfiletype="Neutral Format" -meshfile=two_cubes.neutral -batchmode

meshconvert:
	nmeshimport --netgen two_cubes.neutral two_cubes.nmesh.h5

run:
	nsim two_cubes.py --clean

postprocess1:
	ncol two_cubes 0 m_Co_0 m_Co_1 m_Co_2 m_Py_0 m_Py_1 m_Py_2 | xmgrace -nxy - -param results.par -hdevice PNG -printfile results.png -hardcopy
	mogrify -trim results.png

postprocess2:

	nmagpp --vtk results two_cubes

	#call create vis1.png, vis2, vis3, vis4 manually from mayavi:
	mayavi vis1.mv
	mayavi vis2.mv
	mayavi vis3.mv
	mayavi vis4.mv

trimpngs:
	mogrify -trim vis1.png
	mogrify -trim vis2.png
	mogrify -trim vis3.png
	mogrify -trim vis4.png


scalepngs:
	python ../tools/scalepng cubes.png 10
	python ../tools/scalepng results.png 15
	python ../tools/scalepng vis1.png 8
	python ../tools/scalepng vis2.png 8
	python ../tools/scalepng vis3.png 8
	python ../tools/scalepng vis4.png 8



all: 
	make netgen
	make meshconvert
	make run
	make postprocess1
	make postprocess2 #manual intervention required
	make trimpngs


