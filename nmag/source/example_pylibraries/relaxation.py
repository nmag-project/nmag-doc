import nmag
from nmag import SI, every

# This basic simulation gives us an idea about the relaxation properties
# of the system.

# We find that it takes about 200 frames to do one full oscillation,
# which corresponds to a frequency of about f=1 GHz.

#create simulation object
sim = nmag.Simulation()

# define magnetic material
Py = nmag.MagMaterial(name = 'Py',
                      Ms = SI(1e6, 'A/m'),
                      llg_damping=0.02,
                      exchange_coupling = SI(13.0e-12, 'J/m'))

# load mesh
sim.load_mesh('rod.nmesh.h5',
              [('rod', Py)],
              unit_length = SI(1e-9, 'm'))

# set initial magnetisation
sim.set_m([1,0,0.1])

# set external field
sim.set_H_ext([0,0,0], SI('A/m'))

sim.relax(save =[('fields', every('time', SI(5e-12,"s")))])
