from nmag import every, at
sim.hysteresis(Hs, save=[('averages', 'fields', at('stage_end')),
                         ('restart', at('stage_end') | every(1000, 'step'))])
