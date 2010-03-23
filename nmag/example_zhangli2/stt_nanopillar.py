import nmag, os, math
from nmag import SI, every, at
from nsim.si_units.si import degrees_per_ns

l = 30.0                        # The nanopillar thickness is 30 nm
hl = l/2                        # hl is half the nanopillar thickness
relaxed_m_file = "relaxed_m.h5" # File containing the relaxed magnetisation
mesh_name = "l030.nmesh.h5"     # Mesh name
mesh_unit = SI(1e-9, "m")       # Unit length for space used by the mesh

def run_simulation(sim_name, initial_m, damping, stopping_dm_dt,
                   j, P=0.0, save=[], do=[], do_demag=True):
  # Define the material
  mat = nmag.MagMaterial(
          name="mat",
          Ms=SI(0.8e6, "A/m"),
          exchange_coupling=SI(13.0e-12, "J/m"),
          llg_damping=damping,
          llg_xi=SI(0.01),
          llg_polarisation=P)

  # Create the simulation object and load the mesh
  sim = nmag.Simulation(sim_name, do_demag=do_demag)
  sim.load_mesh(mesh_name, [("np", mat)], unit_length=mesh_unit)

  # Set the pinning at the top and at the bottom of the nanopillar
  def pinning(p):
    x, y, z = p
    tmp = float(SI(x, "m")/(mesh_unit*hl))
    if abs(tmp) >= 0.999:
      return 0.0
    else:
      return 1.0
  sim.set_pinning(pinning)

  if type(initial_m) == str:            # Set the initial magnetisation
    sim.load_m_from_h5file(initial_m)   # a) from file if a string is provided
  else:
    sim.set_m(initial_m)                # b) from function/vector, otherwise

  if j != 0.0:                          # Set the current, if needed
    sim.set_current_density([j, 0.0, 0.0], unit=SI("A/m^2"))

  # Set additional parameters for the time-integration and run the simulation
  sim.set_params(stopping_dm_dt=stopping_dm_dt,
                 ts_rel_tol=1e-7, ts_abs_tol=1e-7)
  sim.relax(save=save, do=do)
  return sim

# If the initial magnetisation has not been calculated and saved into
# the file relaxed_m_file, then do it now!
if not os.path.exists(relaxed_m_file):
  # Initial direction for the magnetisation
  def m0(p):
      x, y, z = p
      tmp = min(1.0, max(-1.0, float(SI(x, "m")/(mesh_unit*hl))))
      angle = 0.5*math.pi*tmp
      return [math.sin(angle), math.cos(angle), 0.0]

  save = [('fields', at('step', 0) | at('stage_end')),
          ('averages', every('time', SI(5e-12, 's')))]

  sim = run_simulation(sim_name="relaxation", initial_m=m0,
                       damping=0.5, j=0.0, save=save,
                       stopping_dm_dt=1.0*degrees_per_ns)
  sim.save_restart_file(relaxed_m_file)
  del sim

# Now we simulate the magnetisation dynamics
save = [('averages', every('time', SI(9e-12, 's')))]
do   = [('exit', at('time', SI(6e-9, 's')))]
run_simulation(sim_name="dynamics", initial_m=relaxed_m_file, damping=0.02,
               j=0.1e12, P=1.0, save=save, do=do, stopping_dm_dt=0.0)
