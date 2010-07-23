# This script is used to compute the equilibrium configuration for the
# magnetisation, which is then used in the second part of the simulation,
# where the dynamics is actually studied.

from thesystem import simulate_nanowire, m0_filename, ps
from nmag.common import *

s = simulate_nanowire(name='relaxation', damping=0.5) # NOTE the high damping!
s.set_m([1, 0, 0])
s.relax(save=[('fields', at('time', 0*ps) | at('convergence'))])
s.save_restart_file(m0_filename)

