netgen:
	netgen -geofile=H.geo -meshfiletype="Neutral Format" -meshfile=H.neutral -batchmode

meshconvert:
	nmeshimport --netgen H.neutral H.nmesh.h5

run1:
	nsim H.py 
	nsim H_novtk.py


post:
	python extract_central_temperature_against_time.py > data_t_maxT.dat

	python extract_T_profile_at_1ns.py > data_x_T_1ns.dat

