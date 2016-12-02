import nmag
from nmag import SI, si, at

# create simulation object
ps = SI(1e-12, "s")

sim = nmag.Simulation()

# define magnetic material (data from Kronmueller's book)
NdFeB = nmag.MagMaterial(name="NdFeB",
                         Ms=1.6*si.Tesla/si.mu0,
                         exchange_coupling=SI(7.3e-12, "J/m"),
                         anisotropy=nmag.uniaxial_anisotropy(axis=[0.01,0.01,1],
                                                             K1=SI(4.3e6, "J/m^3"),
                                                             K2=SI(0*0.65e6, "J/m^3")
                                                             )
                         )

sim.load_mesh("cube.nmesh.h5", [("cube", NdFeB)], unit_length=SI(1.0e-9,"m") )

sim.set_m([-0.01,-0.01,1])

Hs = nmag.vector_set( direction=[0.,0.,1.],
                      norm_list=[-1,-2, [], -4, -4.2, [], -4.9, -4.91, [], -4.94],
                      units=1e6*SI('A/m')
                    )

sim.hysteresis(Hs,
               save=[('fields', 'restart', at('convergence'))]
               )
