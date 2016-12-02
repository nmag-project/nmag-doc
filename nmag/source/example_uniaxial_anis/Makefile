run: simulate media

simulate:
	../../../src/bin/nsim bar.py --clean

media:
	#gnuplot bar_mag_x.gnu
	#gnuplot bar_mag_x_compared.gnu
	python make_plots.py
	python ../tools/scalepng bar_mag_x.png 15
	python ../tools/scalepng bar_mag_x_compared.png 15

