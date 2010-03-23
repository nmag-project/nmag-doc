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
	LHS="\frac{\partial \vec{M}}{\partial t}"; \
	RHS_PREC="-\gamma \, \vec{M} \times \vec{H}"; \
	RHS_DAMP="+\frac{\alpha}{M_{\mathrm{sat}}} \, \vec{M} \times \frac{\partial \vec{M}}{\partial t}"; \
	GRADJ="\hat{j} \cdot \nabla \vec{M}"; \
	RHS_ADIAB="-\frac{v}{M_{\mathrm{sat}}^2} \, \vec{M} \times \left(\vec{M} \times $$GRADJ\right)"; \
	RHS_NADIAB="-\frac{\xi v}{M_{\mathrm{sat}}} \, \vec{M} \times $$GRADJ"; \
	EQ="$$LHS = $$RHS_PREC $$RHS_DAMP \\ $$RHS_ADIAB $$RHS_NADIAB"; \
	$(LATEXPNG) eq_zhangli_llg.png "$$EQ" $(BORDER)
	python ../tools/scalepng eq_zhangli_llg.png 15
	N="\nonumber\\\\"; \
	EQ="\vec{M} & \mathtt{Simulation.set\_m} & \mathrm{is\ the\ magnetisation,} $$N"; \
	EQ="$$EQ \vec{H} & \mathtt{Simulation.set\_H\_ext} & \mathrm{is\ the\ effective\ magnetic\ field,} $$N"; \
	EQ="$$EQ M_{\mathrm{sat}} & \mathtt{Ms} & \mathrm{is\ the\ saturation\ magnetisation,} $$N"; \
	EQ="$$EQ \gamma & \mathtt{llg\_gamma\_G} & \mathrm{is\ the\ gyromagnetic\ ratio,} $$N"; \
	EQ="$$EQ \alpha & \mathtt{llg\_damping} & \mathrm{is\ the\ damping\ parameter,} $$N"; \
	EQ="$$EQ \hat{j} \cdot \nabla &  & \mathrm{derivative\ along\ }\hat{j}"; \
	EQ="$$EQ  \mathrm{,\ the\ direction\ of\ the\ current,} $$N"; \
	$(LATEXPNG) eq_zhangli_llg_d.png "$$EQ" $(BORDER)
	python ../tools/scalepng eq_zhangli_llg_d.png 11
	EQ="v=\frac{P j \mu_B}{e M_{\mathrm{sat}} (1 + \xi^2)}"; \
	$(LATEXPNG) eq_v.png "$$EQ" $(BORDER)
	python ../tools/scalepng eq_v.png 4
	N="\nonumber\\\\"; \
	EQ="j & \mathtt{Simulation.set\_current\_density} & \mathrm{is\ the\ norm\ of\ the\ current\ density,} $$N"; \
	EQ="$$EQ P & \mathtt{llg\_polarisation} & \mathrm{is\ the\ degree\ of\ polarization\ of\ the\ spin\ current,} $$N"; \
	EQ="$$EQ \xi &  \mathtt{llg\_xi} & \mathrm{is\ the\ degree\ of\ non\ adiabadicity.} $$N"; \
	EQ="$$EQ \mu_B &  & \mathrm{is\ the\ Bohr\ magneton,} $$N"; \
	EQ="$$EQ e &  & \mathrm{is\ the\ electron\ charge,} $$N"; \
	$(LATEXPNG) eq_v_d.png "$$EQ" $(BORDER)
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
