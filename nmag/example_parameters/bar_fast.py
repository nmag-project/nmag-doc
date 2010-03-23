import time #need this to measure execution time 
import nmag
from nmag import SI, every, at

mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(0.86e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m"),
                          llg_damping=SI(0.5,""))

sim = nmag.Simulation()

sim.load_mesh("bar30_30_100.nmesh.h5",
              [("Py", mat_Py)],
              unit_length=SI(1e-9,"m"))

sim.set_m([1,0,1])


#Need to call advance time to create timestepper

ps = SI(1e-12, "s") 

sim.advance_time(0*ps)





sim.get_timestepper().set_params(rel_tolerance=1e-4,abs_tolerance=1e-4)

starttime = time.time()

sim.relax(save = [('averages', every('time', 5*ps)),
                  ('fields', at('convergence'))])


stoptime = time.time()

print "Time spent is %5.2f seconds" % (stoptime-starttime)


