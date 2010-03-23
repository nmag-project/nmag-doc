netgen:
	netgen -geofile=block.geo -meshfiletype="Neutral Format" -meshfile=block.neutral -batchmode

meshconvert:
	nmeshimport --netgen block.neutral block.nmesh.h5

run1:
	nsim block.py 


