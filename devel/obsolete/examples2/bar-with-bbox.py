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

sim = nmag.Simulation("bar",mesh_unit_length=SI(1e-9,"m"), fem_only=False )

meshfile = "/tmp/context-problem/bar-with-bbox.nmesh.h5"
sim.load_mesh(meshfile, [("void", []),("Py", mat_Py)])


def initial_magnetization((x, y, z), mag_type):
    import math
    angle_deg = 45
    angle_rad = angle_deg/360.*2*math.pi
    return [math.cos(angle_rad), 0, math.sin(angle_rad)]

sim.set_magnetization(initial_magnetization)



#Once we have a 'restartfield.nvf' file, we can use this to set the initial magnetisation
#sim.set_magnetization('restartfield.nvf')

dt = SI(5e-12, "s") 

target_time = sim.advance_time(intensive_param_by_name, 1e-6*dt)

sim.save_field('m','bar.nvf')

sim.save_data_table()

sim.save_fields_vtk('bar_vtk%04d.vtk' % 0)

for i in range(1, 20):
    time = dt*i
    target_time = sim.advance_time(intensive_param_by_name, time)



    sim.save_fields(fieldnames=sim.fields.keys())

    sim.save_field('m',filename='bar.nvf')

    #write SI data into default file

    sim.save_data_table()
    #sim.save_fields_vtk('sphere-r10_vtk%04d.vtk' % i)

sim.save_fields_vtk('bar_vktfinal.vtk')



