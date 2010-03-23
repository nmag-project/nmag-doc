run: simulate media

simulate:
	../../../../bin/nsim bar.py --clean

media:
	gnuplot bar_mag_x.gnu
	python ../tools/scalepng bar_mag_x.png 15
	gnuplot bar_mag_x_compared.gnu
	python ../tools/scalepng bar_mag_x_compared.png 15

