import nmag
import time
from nmag import SI

# When creating the simulation object, specify that the BEM hmatrix should be
# set up by using the default parameters.
sim = nmag.Simulation(phi_BEM=nmag.default_hmatrix_setup)

# Specify magnetic material, parameters chosen as in example 1
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6, "A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"))

# Load the mesh
sim.load_mesh('sphere.nmesh.h5',
              [('sphere', Py)],
              unit_length=SI(1e-9, 'm'))

# Set the initial magnetisation
sim.set_m([1,0,0])

# Save the demagnetisation field
sim.save_data(fields=['H_demag'])

# Probe the demagnetisation field at ten points within the sphere
for i in range(-5,6):
    x = i*1e-9
    Hdemag = sim.probe_subfield_siv('H_demag', [x,0,0])
    print "x=", x, ": H_demag = ", Hdemag

