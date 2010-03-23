set term postscript eps enhanced color
set out 'hysteresis.eps'
set xlabel 'Applied field H_x 	(A/m)'
set ylabel 'M_x / M_s'
set xrange [-1050000:1050000]
set yrange [-1.2:1.2]
plot 'plot.dat' u 1:2 ti 'ellipsoid example' w lp 3