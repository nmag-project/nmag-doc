.PHONY: default clean

NSIM_ROOT=../../../../..
NSIM_BIN=$(NSIM_ROOT)/bin
NSIM=$(NSIM_BIN)/nsim
NMESHIMPORT=$(NSIM) $(NSIM_BIN)/nmeshimport

nsimfiles_full := $(wildcard thinfilm*_full.py)
nsimfiles_hlib := $(wildcard thinfilm*_hlib.py)
geofiles := $(wildcard thinfilm*.geo)
meshfiles := $(patsubst %.geo, %.nmesh.h5, $(geofiles))

default: ../performance.png ../relaxation_m.png

clean:
	rm -f film*_full_* film*_hlib_* *.dat *.eps *.nmesh.h5 *.neutral \
	*_log.log* convergence.txt test.out *.ps ng.ini

#creating graph comparing the performances with and without hlib

../performance.png: performance.eps
	convert -trim -antialias -density 150 performance.eps \
	-bordercolor white -border 0x0 white performance.png ; \
	python ../../tools/scalepng performance.png 15 ; \
	mv performance.png ../

performance.eps: timings_full.dat timings_hlib.dat
	gnuplot plot_performance_data.gnp

timings_full.dat: $(nsimfiles_full) $(meshfiles)
	./run_simulations_full.sh

timings_hlib.dat: $(nsimfiles_hlib) $(meshfiles)
	./run_simulations_hlib.sh

#create graph comparing m(t) for computation with and without hlib

../relaxation_m.png: relaxation_m.eps
	convert -trim -antialias -density 150 relaxation_m.eps \
	-bordercolor white -border 0x0 relaxation_m.png ; \
	python ../../tools/scalepng relaxation_m.png 15 ; \
	mv relaxation_m.png ../

relaxation_m.eps: relaxation_m.gnp full100.dat hlib100.dat
	gnuplot relaxation_m.gnp

full100.dat: film100_full_dat.ndt
	ncol film100_full_dat.ndt time m_Ni_0 m_Ni_1 m_Ni_2 > full100.dat

film100_full_dat.ndt:
	../../../../../bin/nsim thinfilm100_full.py --clean

hlib100.dat: film100_hlib_dat.ndt
	ncol film100_hlib_dat.ndt time m_Ni_0 m_Ni_1 m_Ni_2 > hlib100.dat

film100_hlib_dat.ndt:
	../../../../../bin/nsim thinfilm100_hlib.py --clean

#creation of meshfile

%.nmesh.h5: %.geo
	netgen -geofile=$< -meshfiletype="Neutral Format" -meshfile=%.neutral \
	-batchmode
	$(NMESHIMPORT) --netgen %.neutral $@ ; \
	gzip -f %.neutral

