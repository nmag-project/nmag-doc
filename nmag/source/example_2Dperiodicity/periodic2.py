import nmag
from nmag import SI, every, at

# define magnetic material
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6,"A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m")
                      )

# size of the simulation cell, plus extra spacing 
x_lattice = 30.01  # the spacing is 0.01 to avoid exchange coupling
y_lattice = 10.01  # between repeated copies of the simualtion cell
z_lattice = 0.0


# list to store the lattice points where the periodic
# copies of the simulation cell will be placed 
lattice_points = []

for xi in range(-4,6):
    for yi in range(-19,21):
        lattice_points.append([xi*x_lattice,yi*y_lattice,0.0*z_lattice])

# create data structure pbc for this macro geometry
pbc = nmag.SetLatticePoints(vectorlist=lattice_points, scalefactor=SI(1e-9,'m'))

#create simulation object, passing macro geometry data structure
sim = nmag.Simulation(periodic_bc=pbc.structure)


# load mesh
sim.load_mesh("prism.nmesh.h5", [("repeated-prism-2D", Py)], unit_length=SI(1e-9,"m") )

# set initial magnetisation
sim.set_m([1.,1.,1.])

# loop over the applied field
s = SI(1,"s")
sim.relax(save=[('averages', 'fields', every('time',5e-12*s) | at('convergence'))])

