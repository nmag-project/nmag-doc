import nmag
from nmag import SI, si

temperature = SI(2.0, "K")
ps = SI(1e-12, "s")
sim = nmag.Simulation(temperature=temperature,
                      user_seed_T = 0,
                      thermal_delta_t=0.001*ps)


# define magnetic material (data from Kronmueller)
NdFeB = nmag.MagMaterial(name="NdFeB",
                         Ms=1.6*si.Tesla/si.mu0,
                         exchange_coupling=SI(7.3e-12, "J/m"),
                         anisotropy=nmag.uniaxial_anisotropy(axis=[0.01,0.01,1],\
                                                             K1=SI(4.3e6, "J/m^3"),\
                                                             K2=SI(0*0.65e6, "J/m^3")))

# load mesh
sim.load_mesh("cube.nmesh.h5", [("cube", NdFeB)], unit_length=SI(1.0e-9,"m") )

# set initial magnetisation from the equilibrium configuration
# with the field = 4.92e6 A/m which has id=19 (check with 'ncol thermal-0K id H_ext_2') 
magn_from_file = nmag.get_subfield_from_h5file('thermal-0K_dat.h5','m_NdFeB',id=10)
sim.set_m(magn_from_file)

# apply external field in -z direction
sim.set_H_ext([0,0,-4.92],unit=SI(1e6,'A/m'))

num_steps = 100
for n in range(0, num_steps):
  print "Step %i/%i" % (n,num_steps)
  sim.advance_time(n*ps)
  if (n % 5) == 0:
    sim.save_data(fields='all')
  else:
    sim.save_data()
