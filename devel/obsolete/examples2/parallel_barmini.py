import nmag
from nmag import SI,mesh
import os

mesh_length=SI(1e-9,"m")

raw_mesh=None

mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m")
                          )

sim = nmag.Simulation(ddd_use_linalg_machine=True)

meshfile = os.path.join("barmini.nmesh.h5")

if os.path.exists(meshfile):
    sim.load_mesh(meshfile, [("PyX", mat_Py)],unit_length=SI(1e-9,"m"))
else:
    sim.defregion("Py", mesh.box([0.0,0.0,0.0], [5, 5, 10]), mat_Py)
    mesh = sim.generate_mesh(([-0.0,-0,-0.0],[5.0,5.0,10.0]), # bounding box
                             a0=2,
                             max_steps=2000,
                             unit_length=mesh_length
                             )
    mesh.save(meshfile)

if ocaml.petsc_mpi_nr_nodes()>1:
    raw_mesh=sim.mesh.raw_mesh
    print "*** RAW MESH *** ",raw_mesh
    ocaml.mesh_set_vertex_distribution(raw_mesh,[100,69])
    sim._init_postmesh()


sim.set_H_ext([0,0,0],SI(1,"A/m"))


def initial_magnetization((x, y, z), mag_type):
    import math
    angle_deg = 45
    angle_rad = angle_deg/360.*2*math.pi
    return [math.cos(angle_rad), 0, math.sin(angle_rad)]

sim.set_magnetization(initial_magnetization)

dt = SI(5e-12, "s") 

target_time = sim.advance_time(1e-6*dt)

sim.save_data_table()

sim.save_fields_vtk('barmini_vtk%04d.vtk' % 0)

target_time = sim.advance_time(5*dt)

sim.save_fields_vtk('barmini_vtk%04d.vtk' % 0)

print sim.get_subfield("m_Py",[2*mesh_length,2*mesh_length,2*mesh_length])
