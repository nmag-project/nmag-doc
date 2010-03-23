# We model a bar 100 nm x 100 nm x 10 nm where a vortex sits in the center.
# This is part I: we just do a relaxation to obtain the shape of the vortex.
import math, nmag
from nmag import SI, at
from nsim.si_units.si import degrees_per_ns

# Define the material
mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(0.8e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m"),
                          llg_gamma_G=SI(0.2211e6, "m/A s"),
                          llg_damping=1.0)

# Define the simulation object and load the mesh
sim = nmag.Simulation()
sim.load_mesh("pyfilm.nmesh.h5", [("Py", mat_Py)], unit_length=SI(1e-9,"m"))

# Set a initial magnetisation which will relax into a vortex
def initial_m(p):
  x, y, z = p
  return [-(y-50.0e-9), (x-50.0e-9), 40.0e-9]

sim.set_m(initial_m)

# Set convergence parameters and run the simulation
sim.set_params(stopping_dm_dt=1.0*degrees_per_ns)
sim.relax(save=[('fields', at('step', 0) | at('stage_end'))])

# Write the final magnetisation to file "vortex_m.h5"
sim.save_restart_file("vortex_m.h5")
