all: updatefiles run media

updatefiles:
	cp ../example1/sphere1.nmesh.h5 .

run:
	../../../../bin/nsim sphere.py --clean

media:
	../../../../bin/nmagpp --vtk=sphere.vtk sphere
	echo "*** Now regenerate .png images from *.mv files! ***"

scalepngs:
	python ../tools/scalepng pinned_core.png 15
	python ../tools/scalepng magnetisation.png 15
