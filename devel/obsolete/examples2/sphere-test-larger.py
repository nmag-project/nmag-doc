import nmag, os, time,sys
from nmag import SI,mesh


mesh_length=SI(1e-9,"m")

clock_time0=time.time()

raw_mesh=None

mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m")
                          )

lam=False

sim = nmag.Simulation(ddd_use_linalg_machine=lam)


meshfile = "/tmp/sphere-test-larger.nmesh.h5"

if os.path.exists(meshfile):
    sim.load_mesh(meshfile, [("sphere", mat_Py)],unit_length=SI(1e-9,"m"))
else:
    sim.defregion("sphere", mesh.ellipsoid([5.0,5.0,8.0]), mat_Py)
    mesh = sim.generate_mesh(([-6.0,-6.0,-10.0],[6.0,6.0,10.0]), # bounding box
                             a0=1.0,
                             max_steps=500,
                             unit_length=mesh_length
                             )
    mesh.save(meshfile)

nr_nodes=ocaml.petsc_mpi_nr_nodes()


if ocaml.petsc_is_mpi():
    print "*** PARALLEL EXECUTION ***"
    raw_mesh=sim.mesh.raw_mesh
    nr_points=ocaml.mesh_nr_points(raw_mesh)
    z=nr_points/nr_nodes
    distrib = [int(round(z*(i+1)))-int(round(z*i)) for i in range(0,nr_nodes)]
    slack=nr_points-reduce(lambda x,y:x+y,distrib)
    distrib[0] = distrib[0] + slack
    print "*** RAW MESH %s *** DISTRIB %s ***" %(repr(raw_mesh),repr(distrib))
    ocaml.mesh_set_vertex_distribution(raw_mesh,distrib)
    sim._init_postmesh()

sim.set_H_ext([0,0,0],SI(1,"A/m"))

def initial_magnetization((x, y, z), mag_type):
    import math
    angle_deg = 45
    angle_rad = angle_deg/360.*2*math.pi
    return [math.cos(angle_rad), 0, math.sin(angle_rad)]

sim.set_magnetization(initial_magnetization)

dt = SI(5e-12, "s") 

clock_time1=time.time()

target_time = sim.advance_time(1e-6*dt)
target_time = sim.advance_time(10*dt)

print "DDD target_time"

result_m=sim.get_subfield("m_Py",[0*mesh_length,0*mesh_length,0*mesh_length])

clock_time2=time.time()

print """Result: %s

Timings (setup/sim/total):  %8.3f sec -- %8.3f sec --  %8.3f sec
"""%(repr(result_m),clock_time1-clock_time0,clock_time2-clock_time1,clock_time2-clock_time0)

sys.stdout.flush()

sim.save_fields_vtk('sphere_test__lam=%s__cpus=%d.vtk'%(repr(lam),nr_nodes))

sys.exit(0)

"""
=== Results ===

data = [0.13872315431435661, 0.4020267364516989, 0.90506073108183105]
data = [0.13881639355586808, 0.40222742896419328, 0.90495729020258686]

Without linalg_machine, old code:

Timings (setup/sim/total):   114.408 sec --  134.183 sec --   248.592 sec

Non-optimized H_demag:

Timings (setup/sim/total):    77.955 sec --  160.679 sec --   238.634 sec

Somewhat optimized H_demag (gmres+ilu, initial_guess_nonzero):

Timings (setup/sim/total):    80.346 sec --  134.815 sec --   215.161 sec


"""
