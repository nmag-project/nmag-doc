# We model a bar 100 nm x 100 nm x 10 nm where a vortex sits in the center.
# This is part II: we load the vortex from file and apply a spin-polarised current

import nmag
from nmag import SI, every, at

# Define the material
mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(0.8e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m"),
                          llg_gamma_G=SI(0.2211e6, "m/A s"),
                          llg_polarisation=1.0,
                          llg_xi=0.05,
                          llg_damping=0.1)

# Define the simulation object and load the mesh
sim = nmag.Simulation()
sim.load_mesh("pyfilm.nmesh.h5", [("Py", mat_Py)], unit_length=SI(1e-9,"m"))

# Set the initial magnetisation: part II uses the one saved by part I
sim.load_m_from_h5file("vortex_m.h5")
sim.set_current_density([1e12, 0, 0], unit=SI("A/m^2"))

sim.set_params(stopping_dm_dt=0.0) # * WE * decide when the simulation should stop!

sim.relax(save=[('fields', at('convergence') | every('time', SI(1.0e-9, "s"))),
                ('averages', every('time', SI(0.05e-9, "s")) | at('stage_end'))],
          do = [('exit', at('time', SI(10e-9, "s")))])
