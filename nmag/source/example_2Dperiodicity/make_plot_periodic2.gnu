set term postscript eps enhanced color 
set out 'periodic2.eps'
set xlabel 'time (s)'
set ylabel 'M / Ms'
plot 'plot_periodic2.dat' u 1:2 ti 'Mx' w lp, 'plot_periodic2.dat' u 1:3 ti 'My' w lp, 'plot_periodic2.dat' u 1:4 ti 'Mz' w lp