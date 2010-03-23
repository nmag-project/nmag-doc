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

sim = nmag.Simulation("sphere", fem_only=False)
unit_length = SI(1e-9,"m")

meshfile = "central-sphere-fine.nmesh.h5"
sim.load_mesh(meshfile, [("Py", mat_Py)], unit_length=unit_length)

def initial_magnetization((x, y, z), mag_type):
    return [1, 0, 0]

sim.set_magnetization(initial_magnetization)

sim.compute_H_fields()
sim.fun_update_energies([])

sim.save_field('m','m-sphere-fem-bem.nvf')
sim.save_data_table()
for p_x in range(-100,101):
    
    try:
        hx,hy,hz = sim.get_H_demag(pos=[p_x*0.1,0.0,0.0],pos_units=unit_length)
        mx,my,mz = sim.get_M([p_x*0.1,0.0,0.0],name="Py",pos_units=unit_length)
    except:
        hx,hy,hz = 0.0,0.0,0.0
        mx,my,mz = 0.0,0.0,0.0
    print "%f %f %f %f %f %f %f" % (p_x*0.1, mx, my,mz, hx, hy, hz)
    
sim.save_fields_vtk('sphere-fem-bem_vtk%04d.vtk' % 0)

