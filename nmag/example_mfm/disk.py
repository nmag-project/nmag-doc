import nmag
from nmag import SI, at

#create simulation object
sim = nmag.Simulation()

# define magnetic material
Py = nmag.MagMaterial( name="Py",
                       Ms=SI(795774,"A/m"),
                       exchange_coupling=SI(13.0e-12, "J/m")
                     )

# load mesh: the mesh dimensions are scaled by 100nm
sim.load_mesh( "../example_mfm/disk.nmesh.h5",
               [("disk", Py)],
               unit_length=SI(1e-9,"m")
             )

# set initial magnetisation
sim.set_m([0.,0.,1.])

# loop over the applied fields Hs
sim.relax(save=[('averages', 'fields', 'restart', at('convergence'))])


