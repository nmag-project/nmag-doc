run:
	make updatefiles
	../../../../bin/nsim sphere_manipulate.py --clean | tee output.txt

vtk:
	nmagpp -v step7.vtk --range 1 sphere_manipulate
	echo "Create step7.png manually from mayavi"

scalepngs:
	python ../tools/scalepng step7.png 15

updatefiles:
	cp ../example1/sphere1.nmesh.h5 .