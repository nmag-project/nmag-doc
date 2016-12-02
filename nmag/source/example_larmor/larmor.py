import nmag
from nmag import SI, every, at, si


#create simulation object and switch off
#the computation of the demagnetising field
sim = nmag.Simulation(do_demag = False)

# define magnetic material so that Js = mu0*Ms = 1 T
Py = nmag.MagMaterial(name="Py",
                      Ms=1.0*si.Tesla/si.mu0,
                      exchange_coupling=SI(13.0e-12, "J/m"),
                      llg_damping = SI(0.0)
                      )
# load mesh
sim.load_mesh("sphere1.nmesh.h5",
              [("sphere", Py)],
              unit_length=SI(1e-9,"m")
             )

# set initial magnetisation
sim.set_m([1,1,1])

# set external field
Hs = nmag.vector_set(direction=[0.,0.,1.],
                     norm_list=[1.0],
                     units=1e6*SI('A/m')
                    )

ps = SI(1e-12, "s")  # ps corresponds to one picosecond

# let the magnetisation precess around the direction of the applied field
sim.hysteresis(Hs,
               save=[('averages', every('time', 0.1*ps))],
               do=[('exit', at('time', 300*ps))])
