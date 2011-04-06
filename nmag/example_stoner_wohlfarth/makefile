
# Note: this does not include "comparisonplot":
all: pngformulae run ncol plot


remesh: netgen meshconvert


#pngformulae:
#	../latex-png demag_nz.png "N_x=\\frac{1}{m^2-1}\\cdot\\left[\\frac{m}{2\\sqrt{m^2-1}}\\cdot\\ln\\left(\\frac{m+\\sqrt{m^2-1}}{m-\\sqrt{m^2-1}}\\right)-1\\right]"
#	python ../tools/scalepng demag_nz.png 11
#	../latex-png Ha.png "H_a=M_s\cdot\Delta N=M_s\cdot\left(N_z-N_x\right)"
#	python ../tools/scalepng Ha.png 7
#	../latex-png hc.png "h_c=\\frac{H_c}{H_a}=\\sin\\theta_0\\cdot\\cos\\theta_0"
#	python ../tools/scalepng hc.png 5


scalepngs:
	python ../tools/scalepng ellipsoid_mesh.png 10


netgen:
	netgen -geofile=ellipsoid.geo -moderate -meshfiletype="Neutral Format" -meshfile=ellipsoid.neu -batchmode

meshconvert:
	nmeshimport --netgen ellipsoid.neu ellipsoid.nmesh.h5
	gzip -f ellipsoid.neu


run:
	../../../../bin/nsim ellipsoid.py --clean

ncol:
	../../../../bin/ncol ellipsoid H_ext_0 H_ext_1 H_ext_2 m_Py_0 m_Py_1 m_Py_2 > plot.dat

plot:   pngformulae
	gnuplot make_plot.gnu
	convert -trim -antialias -density 150 hysteresis.eps hysteresis.png	
	python ../tools/scalepng hysteresis.png 15

clean:
	rm ellipsoid_dat* 

pdf:
	rst2latex --strict --stylesheet=stylesheet.tex -v doc.txt doc.tex
	pdflatex doc.tex; 

