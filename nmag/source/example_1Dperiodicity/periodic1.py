import nmag
from nmag import SI

# define magnetic material
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6,"A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m")
                      )

# size of simulation cell, plus extra spacing
# to avoid exchange interaction across interfaces
# between repeated copies of the simulation cell.
x_lattice = 15.01  # the spacing is 0.01
y_lattice = 0.0
z_lattice = 0.0


# list to store the lattice points where the periodic
# copies will be placed 
lattice_points = []

for xi in range(-3,4):
    lattice_points.append([xi*x_lattice,0.0*y_lattice,0.0*z_lattice])

#  create data structure pbc for this macro geometry
pbc = nmag.SetLatticePoints(vectorlist=lattice_points, scalefactor=SI(1e-9,'m'))

#create simulation object,  passing macro geometry data structure
sim = nmag.Simulation(periodic_bc=pbc.structure)

# load mesh
sim.load_mesh("cube1.nmesh.h5", [("repeated-cube-1D", Py)], unit_length=SI(1e-9,"m") )

# set initial magnetisation along the periodic axis
sim.set_m([1.0,0,0])

# compute the demagnetising field
sim.advance_time(SI(0,"s"))

# probe demag field at the centre of the cube, function
# returns an SI-Value ('siv')
H_demag = sim.probe_subfield_siv('H_demag', [0,0,0])

print "H_demag_x at centre of cube = ", H_demag[0]
print "H_demag_y at centre of cube = ", H_demag[1]
print "H_demag_z at centre of cube = ", H_demag[2]

