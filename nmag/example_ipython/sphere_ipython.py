import nmag
from nmag import SI

# Create simulation object
sim = nmag.Simulation()

# Define magnetic material
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6,"A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"))

# Load mesh
sim.load_mesh("sphere1.nmesh.h5", [("sphere", Py)], unit_length=SI(1e-9,"m"))

# Set initial magnetisation
sim.set_m([1, 0, 0])

# Activate interactive python session
nmag.ipython()

print "Back in main code"
