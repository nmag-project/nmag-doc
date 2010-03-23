import nmag
from nmag import SI, every, at
from numpy import array
import math

# Create simulation object (no demag field!)
sim = nmag.Simulation(do_demag=False)

# Define magnetic material (data from OOMMF materials file)
Co = nmag.MagMaterial(name="Co",
                      Ms=SI(1400e3, "A/m"),
                      exchange_coupling=SI(30e-12, "J/m"),
                      anisotropy=nmag.uniaxial_anisotropy(axis=[0, 0, 1], K1=SI(520e3, "J/m^3")))

# Load the mesh
sim.load_mesh("bar.nmesh.h5", [("bar", Co)], unit_length=SI(1e-9,"m") )

# Our bar is subdivided into 3 regions:
# - region A: for x < offset;
# - region B: for x between offset and offset+length
# - region C: for x > offset+length;
# The magnetisation is defined over all the three regions,
# but is pinned in region A and C.
offset = 2e-9   # m (meters)
length = 500e-9 # m

# Set initial magnetisation
def sample_m0((x, y, z)):
  # relative_position goes linearly from -1 to +1 in region B
  relative_position = -2*(x - offset)/length + 1
  mz = min(1.0, max(-1.0, relative_position))
  return [0, math.sqrt(1 - mz*mz), mz]

sim.set_m(sample_m0)

# Pin magnetisation outside region B
def sample_pinning((x, y, z)):
  return x >= offset and x <= offset + length

sim.set_pinning(sample_pinning)

# Save the magnetisation along the x-axis
def save_magnetisation_along_x(sim):
  f = open('bar_mag_x.dat', 'w')
  for i in range(0, 504):
    x = array([i+0.5, 0.5, 0.5]) * 1e-9
    M = sim.probe_subfield_siv('M_Co', x)
    print >>f, x[0], M[0], M[1], M[2]

# Relax the system
sim.relax(save=[(save_magnetisation_along_x, at('convergence'))])
