import numpy
cut=numpy.fromfile('cut.dat')
xs=numpy.fromfile('xs.dat')

lines = len(cut)
cut2 = cut.reshape(lines/100,100)

assert len(cut2[10]) == len(xs),"internal error (1)"

#true data
#for x,T in map(None,xs,cut2[10]):
#    print "%g\t%g" % (x,T)


for i in range(len(xs)/2):
    print "%g\t%g" % (xs[i],cut2[10][i])

#print "--"
for i in range(len(xs)/2,len(xs)-1):
    print "%g\t%g" % (xs[i],cut2[10][-i-1])


