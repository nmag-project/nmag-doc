import nmag
from nmag import SI,mesh
import nmesh 
import os

H_x = SI(0.0e6, "A/m") 
H_y = SI(0.0e6, "A/m") 
H_z = SI(0.0e6, "A/m") 

intensive_param_by_name={"H_x":H_x, "H_y":H_y, "H_z":H_z}

mat_void = nmag.MagMaterial(name="void",
                          Ms=SI(0.1e6,"A/m"),
                          exchange_coupling=SI(0.0e-12, "J/m")
                          )
mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m")
                          )

sim = nmag.Simulation("sphere",mesh_unit_length=SI(1e-9,"m"), fem_only=False )

meshfile = "sphere-fem-bem.nmesh.h5"
sim.load_mesh(meshfile, [("Py", mat_Py)])

sim.save_mesh("/tmp/debug.mesh")

def initial_magnetization(coords,mag_type):
    if coords[0] > 0:
        return [1, 0, 0]
    else:
        return [-1, 0, 0]

sim.set_magnetization(initial_magnetization)
dt = SI(1.0e-12, "s") 

#sim.advance_time(intensive_param_by_name,SI(1e-16,"s"))
sim.recompute_H(intensive_param_by_name)

#sim.recompute_Energies(intensive_param_by_name)

sim.save_field('m','sphere-fem-bem.nvf')
sim.save_fields(fieldnames=sim.fields.keys())
sim.save_data_table()

# XXX NOTE: get_M seems to be broken!
print "*** M ***",sim.get_m([0.0,0.0,0.0])
print "*** H_d ***",sim.get_H_demag([0.0,0.0,0.0])

print "================================================"
print ocaml.mumag_vector_field_max_angles(sim.fields["m"])
print ocaml.mumag_vector_field_max_angles(sim.fields["H_demag"])
print "================================================"


sim.save_fields_vtk('sphere-fem-bem-vtk%04d.vtk' % 0)
