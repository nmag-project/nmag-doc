import nmag
from nmag import SI

import time #python standard modules, used to measure run time

def run_sim(tol):
    """Function that is called repeatedly with different tolerance values.
    Each function call is carrying out one simulation.
    """
    mat_Py = nmag.MagMaterial(name="Py",
                              Ms=SI(0.86e6,"A/m"),
                              exchange_coupling=SI(13.0e-12, "J/m"),
                              llg_damping=0.5)

    #compose name of simulation to inlude value of tolerance
    sim = nmag.Simulation("bar_%.6f" % tol)

    sim.load_mesh("bar30_30_100.nmesh.h5",
                  [("Py", mat_Py)],
                  unit_length=SI(1e-9,"m"))

    sim.set_m([1,0,1])

    #set tolerance (has to be called after set_m())
    sim.set_params(ts_abs_tol=tol, ts_rel_tol=tol)

    dt = SI(2.5e-12, "s") 

    timing = 0 #initialise variable to measure execution time

    for i in range(0, 121): 
        timing -= time.time()    #start measuring time
        sim.advance_time(dt*i)   #compute time development for 300ps
        timing += time.time()    #stop measuring time
                                 #we exclude time required to save data
                                 
        sim.save_data()          #save averages every 2.5 ps

    #at end of simulation, write performance data into summary file
    f=open('resultsummary.txt','a') #open file to append 
    f.write('%g %d %g\n' % (tol,sim.clock['step'],timing))
    f.close()


#main program
tols = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]

for tol in tols:
    run_sim(tol)


