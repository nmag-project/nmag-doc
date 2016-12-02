import nmag.nmag_multiferroic as nmag
from nmag import SI, every, at

mat_Py = nmag.MagMaterial( name="Py",
                           Ms=SI(0.86e6, "A/m"),
                           exchange_coupling=SI(13.0e-12, "J/m"),
                           multiferroic_coupling = SI(1.0e3,'1/A s'),
                           llg_damping=0.5)


sim = nmag.Simulation("bar_relax")

sim.load_mesh("bar30_30_100.nmesh.h5", [("Py", mat_Py)],
              unit_length=SI(1e-9, "m"))

sim.set_m([1, 0, 1])
sim.set_Electric_ext([1e3,0,0],SI("V/m"))

ps = SI(1e-12,"s")
sim.relax(save = [('averages', every('time', 5*ps)),
                  ('fields', at('convergence'))])

