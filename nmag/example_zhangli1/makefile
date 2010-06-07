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

.PHONY: all clean vtk run pngformulae

all: pngformulae

$(EXAMESH):
	(cd $(CUBICMESH_DIR) && make)

pyfilm.nmesh.h5: $(EXAMESH)
	$(EXAMESH) "pyfilm.nmesh,100:30,100:35,10:3"
	$(NMESHPP) -c pyfilm.nmesh $@
	rm -f pyfilm.nmesh

remesh: netgen meshconvert

pngformulae:
	$(LATEXPNG) -f eq_zhangli_llg.tex eq_zhangli_llg.png $(BORDER)
	python ../tools/scalepng eq_zhangli_llg.png 15
	$(LATEXPNG) -f quant_table1.tex eq_zhangli_llg_d.png $(BORDER)
	python ../tools/scalepng eq_zhangli_llg_d.png 11
	EQ="v=\frac{P j \mu_B}{e M_{\mathrm{sat}} (1 + \xi^2)}"; \
	$(LATEXPNG) eq_v.png "$$EQ" $(BORDER)
	python ../tools/scalepng eq_v.png 4
	$(LATEXPNG) -f quant_table2.tex eq_v_d.png $(BORDER)
	python ../tools/scalepng eq_v_d.png 11

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

clean:
	rm -f *_dat* *_log* *.log *.conf *_progress.txt vortex_m.h5

pdf:
	rst2latex --strict --stylesheet=stylesheet.tex -v doc.txt doc.tex
	pdflatex doc.tex; 
