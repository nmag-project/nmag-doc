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
sim.load_mesh( "nanodot1.nmesh.h5",
               [("cylinder", Py)],
               unit_length=SI(100e-9,"m")
             )

# set initial magnetisation
sim.set_m([1.,0.,0.])

Hs = nmag.vector_set( direction=[1.,0.,0.],
                      norm_list=[1000.0, 900.0, [],
                                 95.0, 90.0, [],
                                 -100.0, -200.0, [],
                                 -1000.0, -900.0, [],
                                 -95.0, -90.0, [],
                                 100.0, 200.0, [], 1000.0],
                      units=1e3*SI('A/m')
                    )


# loop over the applied fields Hs
sim.hysteresis(Hs,
               save=[('averages', 'fields', 'restart', at('convergence'))]
               )

