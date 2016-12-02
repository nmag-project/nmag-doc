set term postscript enhanced color
set out 'no_periodic.ps'
set xlabel 'time (s)'
set ylabel 'M / Ms'
plot 'plot_no_periodic.dat' u 1:2 ti 'Mx' w lp, \
     'plot_no_periodic.dat' u 1:3 ti 'My' w lp, \
     'plot_no_periodic.dat' u 1:4 ti $
