set term png giant size 800, 600
set out 'bar_mag_x.png'
set xlabel 'x (nm)'
set ylabel 'M.z (millions of A/m)'

plot [0:504] [-1.5:1.5] \
  1.4 t "" w l 0, -1.4 t "" w l 0, \
  'bar_mag_x.dat' u ($1/1e-9):($4/1e6) t 'nmag' w l 2
