set term postscript eps enhanced color 
set out 'periodic1_along_periodic_axis.eps'
set xlabel 'copies'
set ylabel 'Demag field'
plot  'demag-field-along-periodic-axis.dat' u 1:3 ti 'nmag' w lp 1 , 'demag-field-along-periodic-axis.dat' u 1:2 ti 'oommf' w p 2

set term postscript eps enhanced color 
set out 'periodic1_out_of_periodic_axis.eps'
set xlabel 'copies'
set ylabel 'Demag field'
plot  'demag-field-out-of-periodic-axis.dat' u 1:3 ti 'nmag' w lp 1 , 'demag-field-out-of-periodic-axis.dat' u 1:2 ti 'oommf' w p 2

