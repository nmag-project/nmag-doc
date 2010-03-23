import nmag
from nmag import SI, every, at

# define magnetic material
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6,"A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m")
                      )

#create simulation object
sim = nmag.Simulation()

# load mesh
sim.load_mesh("prism.nmesh.h5", [("no-periodic", Py)], unit_length=SI(1e-9,"m") )

# set initial magnetisation
sim.set_m([1.,1.,1.])

# loop over the applied field
s = SI(1,"s")

sim.relax(save=[('averages','fields', every('time',5e-12*s) | at('convergence'))])

