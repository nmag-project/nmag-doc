import nmag, os, time,sys
from nmag import SI,mesh


mesh_length=SI(1e-9,"m")

clock_time0=time.time()

raw_mesh=None

mat_Py = nmag.MagMaterial(name="Py",
                          Ms=SI(1e6,"A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m")
                          )

lam=True

sim = nmag.Simulation(ddd_use_linalg_machine=lam)


meshfile = "/tmp/sphere-test.nmesh.h5"

if os.path.exists(meshfile):
    sim.load_mesh(meshfile, [("sphere", mat_Py)],unit_length=SI(1e-9,"m"))
else:
    sim.defregion("sphere", mesh.ellipsoid([5.0,5.0,8.0]), mat_Py)
    mesh = sim.generate_mesh(([-6.0,-6.0,-10.0],[6.0,6.0,10.0]), # bounding box
                             a0=2.0,
                             max_steps=1000,
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

[### ddd_use_linalg_machine=False, 1-cpu ###]

Result: [SI(0.14875680323746421,[]), SI(0.40265109870680693,[]), SI(0.9031871619748596,[])]

Setup Time:     9.830 sec
Sim Time:      10.292 sec

Total Time:    20.122 sec

[### ddd_use_linalg_machine=True, 1-cpu ###]

Result: [SI(0.14876210414171767,[]), SI(0.40265071133584829,[]), SI(0.90318706187280628,[])]

Setup Time:     8.275 sec
Sim Time:      10.049 sec

Total Time:    18.324 sec

[### ddd_use_linalg_machine=True, mpirun -np 2 ###]

Result: [SI(0.14873974881642205,[]), SI(0.40263905632334751,[]), SI(0.90319585298470451,[])]

Setup Time:     9.973 sec
Sim Time:      14.620 sec

Total Time:    24.593 sec

[### ddd_use_linalg_machine=True, mpirun -np 1 ###]

Result: [SI(0.14876210414171767,[]), SI(0.40265071133584829,[]), SI(0.90318706187280628,[])]

Setup Time:    10.282 sec
Sim Time:      10.679 sec

Total Time:    20.961 sec

[### ddd_use_linalg_machine=True, mpirun -np 4 ###]

Hangs after

DDD post make_timestepper

=== TIMINGS ===

Timings (setup/sim/total):     8.014 sec --   10.045 sec --    18.059 sec
Timings (setup/sim/total):     8.336 sec --   10.149 sec --    18.485 sec
Timings (setup/sim/total):     8.249 sec --   10.013 sec --    18.262 sec

After setting KSP's initial_guess_nonzero:

Timings (setup/sim/total):     8.970 sec --   10.345 sec --    19.315 sec
Timings (setup/sim/total):     7.974 sec --   10.373 sec --    18.347 sec

After setting KSP's initial_guess_nonzero, with ksp_type="cg" and pc_type="ilu":

Timings (setup/sim/total):     8.723 sec --    3.919 sec --    12.642 sec

...but completely bogus results!


~~~~~~~~~~~~~

-laplace_phi1, cg, no ILU:

Timings (setup/sim/total):    10.029 sec --   10.290 sec --    20.318 sec
Timings (setup/sim/total):     8.660 sec --   10.154 sec --    18.814 sec

-laplace_phi1, cg, no ILU:

...takes forever!


"""
