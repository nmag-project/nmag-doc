NSIMROOT=../../../..

CUBICMESH_DIR=$(NSIMROOT)/utils/cubicmesh
CUBICMESH=$(CUBICMESH_DIR)/cubicmesh
EXAMESH=$(CUBICMESH_DIR)/examesh
NMESHPP=$(NSIMROOT)/bin/nmeshpp
NMAGPP=$(NSIMROOT)/bin/nmagpp
NMESHIMPORT=$(NSIMROOT)/bin/nmeshimport
NSIM=$(NSIMROOT)/bin/nsim
NCOL=$(NSIMROOT)/bin/ncol
LATEXPNG=python ../latex-png.py
BORDER=20

.PHONY: all clean vtk run mesh

all: run m_of_t.png

$(EXAMESH):
	(cd $(CUBICMESH_DIR) && make)

mesh: l030.nmesh.h5

l030.nmesh.h5: l030.geo
	$(NETGEN) -geofile=l030.geo -moderate -meshfiletype="Neutral Format" -meshfile=l030.neu -batchmode
	$(NMESHIMPORT) --netgen l030.neu l030.nmesh.h5
	rm -f l030.neu

run: l030.nmesh.h5
	$(NSIM) stt_nanopillar.py --clean

dynamics_dat.ndt: run


vtk: relaxation_dat.h5
	rm -f m-*.vtk
	$(NMAGPP) relaxation --vtk=m.vtk

m_of_t.dat: dynamics_dat.ndt
	$(NCOL) dynamics time M_mat_0 M_mat_1 M_mat_2 > $@

m_of_t.eps: m_of_t.gnp m_of_t.dat
	gnuplot m_of_t.gnp

m_of_t.png: m_of_t.eps
	convert -trim -antialias -density 150 m_of_t.eps m_of_t.png
	python ../tools/scalepng m_of_t.png 15

clean:
	rm -f dynamics_*.* relaxation_*.* relaxed_m.h5 *.conf *.log

pdf:
	rst2latex --strict --stylesheet=stylesheet.tex -v doc.txt doc.tex
	pdflatex doc.tex

scale-mesh-png:
	python ../tools/scalepng mesh.png 8
