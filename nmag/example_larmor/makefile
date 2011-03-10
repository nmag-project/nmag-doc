all:
	make pngformulae
	make run
	make ncol
	make plot


updatefiles:
	cp ../example1/sphere1.nmesh.h5 .
run:
	make updatefiles
	../../../../bin/nsim larmor.py --clean
ncol:
	../../../../bin/ncol larmor time m_Py_0 m_Py_1 m_Py_2 > data.txt

plot:   
	gnuplot make_plot.gnu
	python ../tools/scalepng larmor_plot.png 15

clean:
	rm larmor_dat* 

#pngformulae:
#	../latex-png torque.png "{\\mathbf T} = \\mu_0 {\\mathbf m} \\times {\\mathbf H}_{\rm ext}"
#	python ../tools/scalepng torque.png 3.5
#	../latex-png frequency.png "f_{\\rm Larmor}=\\frac{\\gamma}{2\\pi}\\cdot\\mu_0\left|{\\mathbf H}_{\rm ext}\right|"
#	python ../tools/scalepng frequency.png 4.5
#	../latex-png gamma.png "\\gamma=\\frac{g\\cdot e}{2m_e}\\approx 1.7588\\times10^{11}\,{\rm T}^{-1}{\rm s}^{-1}"
#	python ../tools/scalepng gamma.png 7

_rst:
	echo "Only for testing: -> example_larmor.html"
	rst2html example_larmor.txt example_larmor.html
