import nmag
from nmag import SI,mesh
import os

mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m")
                          )

#sim = nmag.Simulation()
sim = nmag.Simulation(ddd_use_linalg_machine=True)

sim.set_H_ext([0,0,0],SI(1,"A/m"))

meshfile = os.path.join("barmini.nmesh.h5")

if os.path.exists(meshfile):
    sim.load_mesh(meshfile, [("PyX", mat_Py)],unit_length=SI(1e-9,"m"))
else:
    sim.defregion("Py", mesh.box([0.0,0.0,0.0], [20, 20, 80]), mat_Py)
    mesh = sim.generate_mesh(([0.0,0.0,0.0],[20.0,20.0,80.0]), # bounding box
                             a0=4.0,
                             max_steps=10,
                             unit_length=SI(1e-9,"m")
                             )
    mesh.save(meshfile)
    

def initial_magnetization(xyz, mag_type):
    import math
    return [math.cos(xyz[2]),math.sin(xyz[2]),0.0]

sim.set_magnetization(initial_magnetization)

for i in range(1,100):
    target_time = sim.advance_time(SI(5e-12*i,"s"),max_it=10000)
    print i, sim.get_H_exchange([5.0,5.0,20.0],pos_unit=SI(1e-9,"m")), sim.get_M([5.0,5.0,20.0],pos_unit=SI(1e-9,"m"))

# With linalg_machine:
#[11453457.390337216, -26135445.98823097, 0.59336900049725783]
#[216782.34319233792, 1108715.0056986893, -0.7466342642118553]

# Without:
#[11453457.390337216, -26135445.98823097, 0.5933690004972646]
#[216782.34319233792, 1108715.0056986893, -0.74663426421185541]
