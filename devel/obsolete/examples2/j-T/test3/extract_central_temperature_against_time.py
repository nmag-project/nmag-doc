import numpy
cut=numpy.fromfile('cut.dat')
xs=numpy.fromfile('xs.dat')

lines = len(cut)
cut2 = cut.reshape(lines/100,100)



for i in range(len(cut2)):
    print "%g\t%g" % (i*0.1e-9,max(cut2[i]))

    
