import time
import nmag
from nmag import SI,mesh
import os, math

msat1 = msat2 = 1e6 * SI("A/m")
mat1 = nmag.MagMaterial(name="M1",
                        Ms=msat1,
                        exchange_coupling=SI(13.0e-12, "J/m"))

mat2 = nmag.MagMaterial(name="M2",
                        Ms=msat2,
                        exchange_coupling=SI(13.0e-12, "J/m"))

sim = nmag.Simulation()

sim.load_mesh("touching.nmesh.h5",
              [("M1", mat1), ("M2", mat2)],
	      unit_length=SI(1e-9,"m"))

def m0(r):
  x, y, z = [ri*1e9 for ri in r]
  return [1, 1, 1]

target_time = sim.advance_time(0)

