import nmag
import time
from nmag import SI

# Create an HMatrix setup object
hms = nmag.HMatrixSetup(nmin=50, eps_aca=1e-5, quadorder=2)

# When creating the simulation object, specify that the BEM hmatrix
# should be set up using the object hms.
sim = nmag.Simulation(phi_BEM=hms)

#specify magnetic material, parameters chosen as in example 1
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6, "A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"))

#load the mesh
sim.load_mesh('sphere.nmesh.h5',
              [('sphere', Py)],
              unit_length=SI(1e-9, 'm'))

#set the initial magnetisation
sim.set_m([1,0,0])

#save the demagnetisation field
sim.save_data(fields=['H_demag'])

#probe the demagnetisation field at ten points within the sphere
for i in range(-5,6):
    x = i*1e-9
    Hdemag = sim.probe_subfield_siv('H_demag', [x,0,0])
    print "x=", x, ": H_demag = ", Hdemag

