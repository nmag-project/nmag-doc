import nmag
from nmag import SI, si

# Create the simulation object
sim = nmag.Simulation()

# Define the magnetic material (data from OOMMF materials file)
Fe = nmag.MagMaterial(name="Fe",
                      Ms=SI(1700e3, "A/m"),
                      exchange_coupling=SI(21e-12, "J/m"),
                      anisotropy=nmag.cubic_anisotropy(axis1=[1, 0, 0],
                                                       axis2=[0, 1, 0],
                                                       K1=SI(48e3, "J/m^3")))

# Load the mesh
sim.load_mesh("cube.nmesh", [("cube", Fe)], unit_length=SI(1e-9, "m"))

# Set the initial magnetisation
sim.set_m([0, 0, 1])

# Launch the hysteresis loop
Hs = nmag.vector_set(direction=[1.0, 0, 0.0001],
                     norm_list=[0, 1, [], 19, 19.1, [], 21, 22, [], 50],
                     units=0.001*si.Tesla/si.mu0)
sim.hysteresis(Hs)
