import nmag
from nmag import SI, si

# Create simulation object
sim = nmag.Simulation()

# Define magnetic material: PermAlloy
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(0.86e6, "A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"))

# Load mesh
sim.load_mesh("sphere1.nmesh.h5", [("sphere", Py)], unit_length=SI(1e-9, "m"))

# Set initial magnetisation to +x direction
sim.set_m([1, 0, 0])

# Pin magnetisation at center in radius of 4e-9m
def pin_at_center((x, y, z)):
  if (x*x + y*y + z*z) < (4e-9)*(4e-9):
    return 0.0 # Inside the 4nm sphere -> pin
  else:
    return 1.0 # Outside -> do not pin

sim.set_pinning(pin_at_center)

# Apply external field in +y direction
unit = 0.5*si.Tesla/si.mu0 # 500mT in A/m
sim.set_H_ext([0*unit, 1*unit, 0*unit])

# Relax the magnetisation
sim.relax()
