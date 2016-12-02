import nmag
from nmag import SI

# Create simulation object
sim = nmag.Simulation()

# Define magnetic material
Py = nmag.MagMaterial(name="Py",
                      Ms=SI(1e6, "A/m"),
                      exchange_coupling=SI(13.0e-12, "J/m"))

# Load mesh
sim.load_mesh("sphere1.nmesh.h5", [("sphere", Py)], unit_length=SI(1e-9, "m"))

# Set initial magnetisation
sim.set_m([1, 0, 0])

# Set external field
sim.set_H_ext([0, 0, 0], SI("A/m"))

# Save and display data in a variety of ways
# Step 1: save all fields spatially resolved

sim.save_data(fields='all')

# Step 2: sample demag field through sphere
for i in range(-10, 11):
    x = i*1e-9                      # position in metres
    H_demag = sim.probe_subfield_siv('H_demag', [x, 0, 0])
    print "x =", x, ": H_demag = ", H_demag

# Step 3: sample exchange field through sphere
for i in range(-10, 11):
    x = i*1e-9                      # position in metres
    H_exch_Py = sim.probe_subfield_siv('H_exch_Py', [x, 0, 0])
    print "x =", x, ": H_exch_Py = ", H_exch_Py


# Now modify the magnetisation at position (0,0,0) (this happens to be
# node 0 in the mesh) in steps 4 to 6:

# Step 4: request a vector with the magnetisation of all sites in the mesh
myM = sim.get_subfield('m_Py')

# Step 5:  We modify the first entry:
myM[0] = [0, 1, 0]

# Step 6: Set the magnetisation to the new (modified) values
sim.set_m(myM)

# Step 7: saving the fields again (so that we can later plot the demag
# and exchange field)
sim.save_data(fields='all') 

# Step 8: sample demag field through sphere (as step 2)
for i in range(-10, 11):
    x = i*1e-9                      # position in metres
    H_demag = sim.probe_subfield_siv('H_demag', [x, 0, 0])
    print "x =", x, ": H_demag = ", H_demag

# Step 9: sample exchange field through sphere (as step 3)
for i in range(-10, 11):
    x = i*1e-9                      # position in metres
    H_exch_Py = sim.probe_subfield_siv('H_exch_Py', [x, 0, 0])
    print "x =", x, ": H_exch_Py = ", H_exch_Py
