import nmag
from nmag import SI, every, at

mat_Py1 = nmag.MagMaterial(name="Py1",
                          Ms=SI(0.86e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m"),
                          llg_damping=0.5)

mat_Py2 = nmag.MagMaterial(name="Py2",
                          Ms=SI(0.86e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m"),
                          llg_damping=0.5)

#Example 1: both materials have same magnetisation

sim = nmag.Simulation()

sim.load_mesh("sphere_in_box.nmesh.h5",
              [("Py1", mat_Py1),("Py2", mat_Py2)],
              unit_length=SI(15e-9,"m"))

sim.set_m([1,0,0])

sim.relax(save = [('averages','fields','restart', every('step', 1000) | at('convergence'))])
