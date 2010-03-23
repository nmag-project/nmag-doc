import nmag
from nmag import SI
import math

# define magnetic material
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6,"A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"),
                      llg_damping = SI(0.02,"")
                      )

# lattice spacings along the main axes; 
# the value must be zero for no periodic copies,
# equal to the mesh dimension along the 
# given axis otherwise
x_lattice = 30.0
y_lattice = 0.0
z_lattice = 0.0


# list to store the lattice points where the periodic
# copies will be placed 
lattice_points = []

for xi in range(-1,2):
    lattice_points.append([xi*x_lattice,0.0*y_lattice,0.0*z_lattice])

# copies of the system along the x-axis
pbc = nmag.SetLatticePoints(vectorlist=lattice_points, scalefactor=SI(1e-9,'m'))

#create simulation object
sim = nmag.Simulation(periodic_bc=pbc.structure)

# load mesh
sim.load_mesh("periodic.nmesh", [("periodic-film", Py)], unit_length=SI(1e-9,"m") )

print ocaml.mesh_plotinfo_periodic_points_indices( sim.mesh.raw_mesh )

# function to set the magnetisation 
def perturbed_magnetisation(pos):
    x,y,z = pos
    newx = x*1e9
    newy = y*1e9
    if 8<newx<14 and -3<newy<3:
        # the magnetisation is twisted a bit 
        return [1.0, 5.*(math.cos(math.pi*((newx-11)/6.)))**3 *\
                        (math.cos(math.pi*(newy/6.)))**3, 0.0]
    else:
        return [1.0, 0.0, 0.0]


# set initial magnetisation
sim.set_m(perturbed_magnetisation)

# let the system relax generating spin waves
s = SI("s")
from nsim.when import every, at
sim.relax(save=[('averages','fields', every('time', 0.05e-12*s) | at('convergence'))],
          do=[('exit', at('time', 10e-12*s))])

