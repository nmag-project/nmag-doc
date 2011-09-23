from nmag import every, at
sim.hysteresis(Hs, save=[('averages', 'fields', at('stage_end')),
                         ('restart', at('stage_end') | every('step', 1000))])
