import numpy
cut=numpy.fromfile('cut.dat')
xs=numpy.fromfile('xs.dat')

lines = len(cut)
cut2 = cut.reshape(lines/100,100)

import pylab



for i in range(0,81,10):
    pylab.plot(xs,cut2[i],label='i=%3d, t=%g' % (i,i*0.1e-9))


pylab.ylabel('T [K]')
pylab.grid()

pylab.savefig('lineplot.eps')
pylab.savefig('lineplot.png')

#pylab.show()
