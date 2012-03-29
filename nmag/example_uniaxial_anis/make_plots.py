
#set term png giant size 800, 600
#set out 'bar_mag_x_compared.png'
#
#set xlabel 'x (nm)'
#set ylabel 'M.z (millions of A/m)'

from numpy import arctan, sinh, cos, pi, sqrt
atan=arctan

def Mz(x):
    return 1400e3 * cos(pi/2 + atan(sinh((x - 252e-9)/sqrt(30e-12/520e3))))

import pylab
import numpy as np
x=np.linspace(220e-9,280e-9,100)
Mz0 = Mz(x)


tmp=np.loadtxt('bar_mag_x.dat')
x2 = tmp[:,0]
Mz2= tmp[:,3]


pylab.plot(x2,Mz2,'-',label='nmag')
pylab.xlabel('x [m]')
pylab.ylabel('Mz [A/m]')
#pylab.axis([min(x2),max(x2),min(Mz2)*1.1,max(Mz2)*1.1])
pylab.legend()
pylab.grid()
pylab.savefig('bar_mag_x.png')
pylab.close()


#read OOMMF data
tmp=np.loadtxt('oommf/bar_mag_x.txt')
x3 = tmp[:,0]
Mz3= tmp[:,3]


pylab.plot(x,Mz0,label='analytical')
pylab.plot(x2,Mz2,'o',label='nmag')
pylab.plot(x2,Mz2,'x',label='oommf')
pylab.xlabel('x [m]')
pylab.ylabel('Mz [A/m]')
pylab.axis([min(x),max(x),min(Mz0)*1.1,max(Mz0)*1.1])
pylab.legend()
pylab.grid()
pylab.savefig('bar_mag_x_compared.png')

pylab.show()

#plot [220:280] [-1.5:1.5] \
#  1.4 t "" w l 0, -1.4 t "" w l 0, \
#  'bar_mag_x.dat' u ($1/1e-9):($4/1e6) t 'nmag' w lp 2, \
#  'oommf/bar_mag_x.txt'u ($1/1e-9):($4/1e6) t 'oommf' w lp 1, \
#  Mz(x*1e-9)/1e6 ti 'analytical' w l 3
