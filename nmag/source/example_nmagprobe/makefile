NSIMROOT=../../../src

CUBICMESH_DIR=$(NSIMROOT)/utils/cubicmesh
CUBICMESH=$(CUBICMESH_DIR)/cubicmesh
EXAMESH=$(CUBICMESH_DIR)/examesh
NMESHPP=$(NSIMROOT)/bin/nmeshpp
NMAGPP=$(NSIMROOT)/bin/nmagpp
NMAGPROBE=$(NSIMROOT)/bin/nmagprobe
NMESHIMPORT=$(NSIMROOT)/bin/nmeshimport
NSIM=$(NSIMROOT)/bin/nsim
NCOL=$(NSIMROOT)/bin/ncol
LATEXPNG=python ../latex-png.py
BORDER=20

.PHONY: all clean vtk run

all: remesh relax


netgen:
	netgen -geofile=cylinder.geo -meshfiletype="Neutral Format" -meshfile=cylinder.neutral -batchmode

meshconvert:
	$(NMESHIMPORT) --netgen cylinder.neutral cylinder.nmesh.h5
	rm cylinder.neutral


remesh: netgen meshconvert

relax:
	$(NSIM) relaxation.py --clean

dynamics:
	$(NSIM) dynamics.py --clean

nmagprobe:
	$(NMAGPROBE) --verbose dynamics_dat.h5 --field=m_Py \
	  --time=0,100e-12,101 --space=-300,300,201/0/0 --ref-time=0.0 \
	  --scalar-mode=component,1 --out=real-space.dat \
	  --ft-axes=0,1 --ft-out=norm --ft-out=rec-space.dat

pngs:
	gnuplot plot.gnp

recreal-space:
	gnuplot plot.gnp
	python ../tools/scalepng real-space.png 9
	python ../tools/scalepng rec-space.png 9

clean:
	rm -f *_dat* *_log* *.log *.conf *_progress.txt vortex_m.h5

