netgen:
	netgen -geofile=H.geo -meshfiletype="Neutral Format" -meshfile=H.neutral -batchmode

meshconvert:
	nmeshimport --netgen H.neutral H.nmesh.h5

run1:
	nsim H.py 


