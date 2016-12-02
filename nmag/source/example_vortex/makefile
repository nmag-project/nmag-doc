all:
	make run_nsim
	make ncol
	make plot

run_nsim:
	../../../../bin/nsim nanodot.py --clean

#run_magpar:
#	magpar

ncol:
	../../../../bin/nsim ../../../../bin/ncol nanodot H_ext_0 m_Py_0 > nmag.dat
	#../../../../bin/nsim ../../../../bin/ncol nanodot.log 0 2 3 6 7 8  | grep -v '#' > magpar.dat 

plot:
	gnuplot make_comparison_plot.gnu
	convert -trim -antialias -density 150 nanodot_evo.eps nanodot_evo.png
	python ../tools/scalepng nanodot_evo.png 15


scalepngs:
	python ../tools/scalepng nanodot-1.png 10
	python ../tools/scalepng nanodot-2.png 10
	python ../tools/scalepng nanodot-3.png 10
	python ../tools/scalepng nanodot-4.png 10
	python ../tools/scalepng nanodot-5.png 10
	python ../tools/scalepng nanodot-6.png 10
