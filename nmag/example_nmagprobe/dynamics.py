from thesystem import simulate_nanowire, m0_filename, ps, nm
from nmag.common import *

# Details about the pulse
pulse_boundary = -300.0e-9 + 0.5e-9 # float in nm
pulse_direction = [0, 1, 0]
pulse_amplitude = SI(1e5, 'A/m')
pulse_duration = 1*ps

# Function which sets the magnetisation to zero
def switch_off_pulse(sim):
  sim.set_H_ext([0.0, 0.0, 0.0], unit=pulse_amplitude)

# Function which sets the pulse as a function of time/space
def switch_on_pulse(sim):
  def H_ext(r):
    if r[0] < pulse_boundary:
      return pulse_direction
    else:
      return [0.0, 0.0, 0.0]

  sim.set_H_ext(H_ext, unit=pulse_amplitude)

# Here we run the simulation: do=[....] is used to set the pulse
#   save=[...] is used to save the data.
s = simulate_nanowire('dynamics', 0.05)
s.load_m_from_h5file(m0_filename)
s.relax(save=[('fields', every('time', 0.5*ps))],
        do=[(switch_on_pulse, at('time', 0*ps)),
            (switch_off_pulse, at('time', pulse_duration)),
            ('exit', at('time', 200*ps))])

