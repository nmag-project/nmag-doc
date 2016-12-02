import nmag
from nmag import SI, at

#create simulation object
sim = nmag.Simulation()

# define magnetic material
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6,"A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"))

# load mesh: the mesh dimensions are scaled by 0.5 nm
sim.load_mesh("ellipsoid.nmesh.h5",
              [("ellipsoid", Py)],
              unit_length=SI(1e-9,"m"))

# set initial magnetisation
sim.set_m([1.,1.,0.])

Hs = nmag.vector_set(direction=[1.,1.,0.],
                     norm_list=[1.0, 0.995, [], -1.0,
		                -0.995, -0.990, [], 1.0],
                     units=1e6*SI('A/m'))

# loop over the applied fields Hs
sim.hysteresis(Hs, save=[('averages', at('convergence'))])

