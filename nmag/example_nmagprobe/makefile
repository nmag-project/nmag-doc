NSIMROOT=../../../..

CUBICMESH_DIR=$(NSIMROOT)/utils/cubicmesh
CUBICMESH=$(CUBICMESH_DIR)/cubicmesh
EXAMESH=$(CUBICMESH_DIR)/examesh
NMESHPP=$(NSIMROOT)/bin/nmeshpp
NMAGPP=$(NSIMROOT)/bin/nmagpp
NSIM=$(NSIMROOT)/bin/nsim
NCOL=$(NSIMROOT)/bin/ncol
LATEXPNG=python ../latex-png.py
BORDER=20

.PHONY: all clean vtk run

all:

pyfilm.nmesh.h5: $(EXAMESH)
	$(EXAMESH) "pyfilm.nmesh,100:30,100:35,10:3"
	$(NMESHPP) -c pyfilm.nmesh $@
	rm -f pyfilm.nmesh

remesh: netgen meshconvert

vortex_m.h5:
	$(NSIM) relaxation.py --clean

run: vortex_m.h5
	$(NSIM) stt.py --clean

vtk: vortex_m.h5
	rm -f m-*.vtk
	$(NMAGPP) relaxation --vtk=m.vtk

m_of_t.dat: stt_dat.ndt
	$(NCOL) stt time M_Py_0 M_Py_1 M_Py_2 > $@

m_of_t.eps: m_of_t.gnp m_of_t.dat
	gnuplot m_of_t.gnp

m_of_t.png: m_of_t.eps
	convert -trim -antialias -density 150 m_of_t.eps m_of_t.png
	python ../tools/scalepng m_of_t.png 15

recreal-space:
	gnuplot plot.gnp
	python ../tools/scalepng real-space.png 9
	python ../tools/scalepng rec-space.png 9

clean:
	rm -f *_dat* *_log* *.log *.conf *_progress.txt vortex_m.h5

pdf:
	rst2latex --strict --stylesheet=stylesheet.tex -v doc.txt doc.tex
	pdflatex doc.tex; 
