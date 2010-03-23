from sys import path
from os import getcwd

cwd = getcwd()
if ((cwd in path)==False):
    path.insert(0,cwd)
import simtools

import nmag

hmatrix = nmag.HMatrixSetup()
stats = simtools.run_simulation(simname="film40_hlib",
                                meshfile="thinfilm40.nmesh.h5",
                                hmatrix=hmatrix,tol=0.00001)

number_nodes = stats[0]
mem_rss = stats[1]
T_setup = stats[2]
T_sim = stats[3]

ofile = open('timings_hlib.dat',"a")
ofile.write("%6d   %6.3f   %8.1f   %8.1f\n" % (number_nodes,mem_rss,T_setup,T_sim))
ofile.close()
