#set term postscript enhanced eps color
set term png giant size 800, 600
set out 'bar_mag_x_compared.png'

set xlabel 'x (nm)'
set ylabel 'M.z (millions of A/m)'

Mz(x) = 1400e3 * cos(pi/2 + atan(sinh((x - 252e-9)/sqrt(30e-12/520e3))))

plot [220:280] [-1.5:1.5] \
  1.4 t "" w l 0, -1.4 t "" w l 0, \
  'bar_mag_x.dat' u ($1/1e-9):($4/1e6) t 'nmag' w lp 2, \
  'oommf/bar_mag_x.txt'u ($1/1e-9):($4/1e6) t 'oommf' w lp 1, \
  Mz(x*1e-9)/1e6 ti 'analytical' w l 3
