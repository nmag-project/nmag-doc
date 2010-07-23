import math
from nmag.common import *
from thesystem import simulate_nanowire, m0_filename, ps, nm

# Details about the disturbance region (a sphere)
disturb_origin = [-50.0, 0.0, 0.0]
disturb_thickness = 0.5
disturb_direction = [0, 1, 0]
disturb_amplitude = SI(1e5, 'A/m')
disturb_duration = 1*ps

# Function which sets the magnetisation to zero
def set_to_zero(sim):
  sim.set_H_ext([0.0, 0.0, 0.0], unit=disturb_amplitude)

# Function which sets the disturbance as a function of time/space
c = float(nm/SI('m')) # conversion factor
def set_disturbance(sim):
  def H_ext(r):
    # x, y and z are (floats) in nanometers
    x, y, z = [xi/c - x0i
               for xi, x0i in zip(r, disturb_origin)]
    if x < disturb_thickness:
        return disturb_direction
    else:
        return [0.0, 0.0, 0.0]

  sim.set_H_ext(H_ext, unit=disturb_amplitude)

# Here we run the simulation: do=[....] is used to set the disturbance
#   save=[...] is used to save the data.
s = simulate_nanowire('dynamics', 0.05)
s.load_m_from_h5file(m0_filename)
s.relax(save=[('fields', every('time', 0.5*ps))],
        do=[(set_disturbance, at('time', 0*ps)),
            (set_to_zero, at('time', disturb_duration)),
            ('exit', at('time', 200*ps))])

