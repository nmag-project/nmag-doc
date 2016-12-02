import nmag
from nmag import SI, every, at
from nsim.si_units import si
import math

# Create simulation object (no demag field!)
sim = nmag.Simulation(do_demag=False)

# Function to compute the scalar product of the vectors a and b
def scalar_product(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

# Here we define a function which returns the energy for a uniaxial
# anisotropy of order 4.
K1 = SI(43e3, "J/m^3")
K2 = SI(21e3, "J/m^3")
axis = [0, 0, 1]        # The (normalised) axis
def my_anisotropy(m):
    a = scalar_product(axis, m)
    return -K1*a**2 - K2*a**4

my_material = nmag.MagMaterial(name="MyMat",
                               Ms=SI(1e6, "A/m"),
                               exchange_coupling=SI(10e-12, "J/m"),
                               anisotropy=my_anisotropy,
                               anisotropy_order=4)

# Load the mesh
sim.load_mesh("coin.nmesh.h5", [("coin", my_material)], unit_length=SI(1e-9, "m"))

# Set the magnetization
sim.set_m([-1, 0, 0])

# Compute the hysteresis loop
Hs = nmag.vector_set(direction=[1.0, 0, 0.0001],
                     norm_list=[-0.4, -0.35, [], 0, 0.005, [], 0.15, 0.2, [], 0.4],
                     units=si.Tesla/si.mu0)

sim.hysteresis(Hs, save=[('fields', 'averages', at('convergence'))])
