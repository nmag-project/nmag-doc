
mesh:
	netgen -geofile=thermal-cube.geo -moderate -meshfiletype="Neutral Format" -meshfile=thermal-cube.neu -batchmode
	../../../../bin/nsim ../../../../bin/nmeshimport --netgen thermal-cube.neu thermal-cube.nmesh
	../../../../bin/nsim ../../../../bin/nmeshpp --scale=0.4 thermal-cube.nmesh cube.nmesh.h5

run0K:
	../../../../bin/nsim thermal-0K.py

run2K:
	../../../../bin/nsim thermal-2K.py

plot:
	../../../../bin/ncol thermal-2K time m_NdFeB_2 > t_Mz.dat
	#plot, for example with 
	#xmgrace  -nxy t_Mz.dat

scalepngs:
	python ../tools/scalepng random-noise-netgen-mesh.png 10	
	python ../tools/scalepng netgen-mesh.png 10
	python ../tools/scalepng switching-field-4_92-nmag-netgen.png 15
	python ../tools/scalepng switch-snapshot-2K-netgen-mesh-4_92.png 10
