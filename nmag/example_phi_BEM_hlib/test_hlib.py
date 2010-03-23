import nmag
import time
from nmag import SI

import os



def run_sim(sim):


    # Specify magnetic material, parameters chosen as in example 1
    Py = nmag.MagMaterial(name="Py",
                          Ms=SI(1e6, "A/m"),
                          exchange_coupling=SI(13.0e-12, "J/m"))
    # Load the mesh
    sim.load_mesh('sphere.nmesh.h5',
                  [('sphere', Py)],
                  unit_length=SI(1e-9, 'm'))
    
    # Set the initial magnetisation
    sim.set_m([1,0,0])
    
    # Save the demagnetisation field
    sim.save_data(fields=['H_demag'])
    
    # Probe the demagnetisation field at ten points within the sphere

    results = []
    
    for i in range(-5,6):
        x = i*1e-9
        Hdemag = sim.probe_subfield_siv('H_demag', [x,0,0])
        results.append(Hdemag)
        print "x=", x, ": H_demag = ", Hdemag

    return results


def test_hlib():

    # When creating the simulation object, specify that the BEM hmatrix should be
    # set up by using the default parameters.

    #work out in which directory the data files are
    org_dir = os.getcwd()
    os.chdir(os.path.split(__file__)[0])

    #delete old data files to avoid failure due to this
    def remove_if_exists(filename):
        if os.path.exists(filename):
            os.remove(filename)

    remove_if_exists("plain_dat.ndt")
    remove_if_exists("plain_dat.h5")
    remove_if_exists("hlib_dat.ndt")
    remove_if_exists("hlib_dat.h5")

    sim_plain = nmag.Simulation(name='plain')
    sim_hlib  = nmag.Simulation(name='hlib',phi_BEM=nmag.default_hmatrix_setup)

    res_plain = run_sim(sim_plain) 
    res_hlib  = run_sim(sim_hlib) 

    #extract x-components of results
    hdx1 = [ h[0] for h in res_plain ]
    hdx2 = [ h[0] for h in res_hlib ]

    #Output as of 13 December 2009, svn  6453, HF 
    #
    #for r1,r2 in zip(hdx1,hdx2):
    #   print r1-r2,(r1-r2)/r1,r1,r2
    #
    # produces:
    #    -1.10303013207 3.31174926024e-06 -333065.714037 -333064.611006
    #    -1.18373466324 3.55403967073e-06 -333067.374849 -333066.191114
    #    -1.26686538779 3.80361435718e-06 -333068.831071 -333067.564206
    #    -1.32740328263 3.98536687005e-06 -333069.28218 -333067.954776
    #    -1.34173143021 4.02840496219e-06 -333067.663951 -333066.32222
    #    -1.35605957796 4.07144347306e-06 -333066.045723 -333064.689663
    #    -1.37001170352 4.11335588419e-06 -333064.228355 -333062.858343
    #    -1.34486746625 4.03791095727e-06 -333060.208727 -333058.863859
    #    -1.30307299155 3.91246719093e-06 -333056.592671 -333055.289598
    #    -1.26131784206 3.78713877786e-06 -333052.976414 -333051.715096
    #    -1.19995578279 3.60291609527e-06 -333051.270434 -333050.070478

    #the actual test
    for r1,r2 in zip(hdx1,hdx2):
        assert  (r1-r2)/float(r1)<1e-5

    #this avoids getting relative paths in the checking progress output
    os.chdir(org_dir)


if __name__ == "__main__":
    test_hlib()
    
