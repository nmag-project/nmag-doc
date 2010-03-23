import nmag
from nmag import SI, every, at

sim = nmag.Simulation()

# define magnetic material Cobalt (data from OOMMF materials file)
Co = nmag.MagMaterial(name="Co",
                      Ms=SI(1400e3, "A/m"),
                      exchange_coupling=SI(30e-12, "J/m"),
                      anisotropy=nmag.uniaxial_anisotropy(axis=[0,0,1], K1=SI(520e3, "J/m^3")))

# define magnetic material Permalley
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(860e3,"A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"))

# load mesh
sim.load_mesh("two_cubes.nmesh.h5",
              [("cube1", Py),("cube2", Co)],
              unit_length=SI(1e-9,"m")
              ) 

# set initial magnetisation along the 
# positive x axis for both cubes, slightly off in z-direction
sim.set_m([0.999847695156, 0, 0.01745240643731])

ns = SI(1e-9, "s") # corresponds to one nanosecond

sim.relax(save = [('averages', every('time', 0.01*ns)),
                  ('fields', every('time', 0.05*ns) | at('convergence'))])

