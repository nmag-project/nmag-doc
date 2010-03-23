.PHONY: default clean run_films rerun

NSIM_BIN=../../../../bin
NSIM=$(NSIM_BIN)/nsim
NMESHIMPORT=$(NSIM) $(NSIM_BIN)/nmeshimport

default: sphere_hlib_dat.h5 sphere_hlib2_dat.h5 run_films 

rerun:
	make clean default

run_films:
	(cd simulation_films && make)

sphere_hlib_dat.h5: sphere_hlib.py sphere.nmesh.h5
	$(NSIM) sphere_hlib.py --clean

sphere_hlib2_dat.h5: sphere_hlib2.py sphere.nmesh.h5
	$(NSIM) sphere_hlib2.py --clean

sphere.nmesh.h5: sphere.geo
	netgen -geofile=$< -meshfiletype="Neutral Format" \
	       -meshfile=sphere.neutral -batchmode
	$(NMESHIMPORT) --netgen sphere.neutral $@ ; \
	  gzip -f sphere.neutral

clean:
	rm -f sphere_hlib_* sphere_hlib2_* ng.ini supermatrix.ps \
	   convergence.txt memory_info.dat sphere.nmesh.h5 sphere.neutral.gz
	(cd simulation_films && make clean)

scale_pngs:
	python ../tools/scalepng relaxation_m.png 15
	python ../tools/scalepng thinfilm130.png 15
