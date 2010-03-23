set term postscript eps enhanced color	
set out 'nanodot_evo.eps'
set xlabel 'Applied field (kA/m)'
set ylabel 'M / M_s'
plot [-250:50] [-1.2:1.2] 'magpar.dat' u 2:3 ti 'magpar' w lp 4 , 'nmag.dat' u ($1/1000):2  ti 'nmag' w lp 3

