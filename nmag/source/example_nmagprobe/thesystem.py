# In this file we define the material parameters and geometry of the system
# so that we can use it in two simulations: first during the relaxation,
# then during the dynamics

from nmag.common import *

nm = SI(1e-9, 'm')    # define nm as "nanometre"
ps = SI(1e-12, 's')   # define ps as "picosecond"

m0_filename = "m0.h5" # the file containing the equilibrium magnetisation

# A function which sets up the simulation with given name and damping
def simulate_nanowire(name=None, damping=0.5):
  permalloy = MagMaterial('Py',
                          Ms=SI(0.86e6, 'A/m'),
                          exchange_coupling=SI(13e-12, 'J/m'),
                          llg_damping=damping)

  s = Simulation(name)
  s.load_mesh("cylinder.nmesh.h5", [('nanopillar', permalloy)],
              unit_length=nm)
  return s

