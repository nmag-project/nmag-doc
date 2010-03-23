import nmag
from nmag import SI,mesh
import nmesh 
import os

H_x = SI(0.0e6, "A/m") 
H_y = SI(0.0e6, "A/m") 
H_z = SI(0.0e6, "A/m") 

intensive_param_by_name={"H_x":H_x, "H_y":H_y, "H_z":H_z}

mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m")
                          )

sim = nmag.Simulation("sphere", fem_only=True)

meshfile = "bar-with-bbox.nmesh.h5"
sim.load_mesh(meshfile, [("void",[]),("Py", mat_Py)],unit_length=SI(1e-9,"m"))


def initial_magnetization((x, y, z), mag_type):
    return [1, 0, 0]

sim.set_magnetization(initial_magnetization)

sim.compute_H_fields()
sim.fun_update_energies([])

sim.save_field('m','m-sphere-fem-only.nvf')
sim.save_data_table()

sim.save_fields_vtk('sphere_vtk%04d.vtk' % 0)

