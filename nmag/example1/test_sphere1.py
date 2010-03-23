"""This test program is called by py.test (from many directory levels above)."""

import os

def do_the_simulation(the_sim=[]):

    #use 'the_sim' to run the simulation only once
    if len(the_sim) == 1:
        return the_sim[0]

    #py.test will pass on any command line arguments to the code it tests.
    #This will confuse nmag (as it doesn't know how to handle py.test's '-k'
    #and '--verbose' switches, etc.
    #
    #We thus manually delete all other entries that the name of the executable
    #from sys.argv.
    import sys

    proper_sysargv = sys.argv
    sys.argv = sys.argv[0:1]

    #work out in which directory the data files are
    org_dir = os.getcwd()
    os.chdir(os.path.split(__file__)[0])

    #delete old data files to avoid failure due to this
    def remove_if_exists(filename):
        if os.path.exists(filename):
            os.remove(filename)

    remove_if_exists("pytest_main_dat.ndt")
    remove_if_exists("pytest_main_dat.h5")

    #import the example from the manual (i.e. execute main program in there)
    import sphere1

    #this avoids getting relative paths in the checking progress output
    os.chdir(org_dir)


    the_sim.append(sphere1)

    return sphere1


def test_demag_at_origin():

    sphere1 = do_the_simulation()

    #probe H_demag at origin, should be Ms/3
    H_demag = sphere1.sim.probe_subfield_siv('H_demag', [0,0,0])

    tol = 0.01

    Ms = sphere1.Py.Ms.value

    H_demag_analytic = -Ms/3.

    assert abs(H_demag[0]/H_demag_analytic-1) < tol, \
           "H_demag is %s, expect %s, tolerance is %s" % \
           (H_demag[0],H_demag_analytic,tol)



def test_get_subfield_average_type():

    sphere1 = do_the_simulation()

    M_average = sphere1.sim.get_subfield_average_siv("M","Py")

    assert type(M_average) == list, "return value of get_subfield_average is not list"


def test_get_subfield_average_value():

    sphere1 = do_the_simulation()

    M_average = sphere1.sim.get_subfield_average_siv("M","Py")

    Ms = sphere1.Py.Ms.value
    assert abs(M_average[0] - Ms)/Ms < 1e-15 #this should be zero to machine accuracy


def test_get_subfield_commands():

    sphere1 = do_the_simulation()

    n_nodes = len(sphere1.sim.mesh.points)
    dummy = sphere1.sim.get_subfield("M_Py")
    assert len(dummy) == n_nodes
    pos = sphere1.sim.get_subfield_positions("M_Py")
    assert len(dummy) == n_nodes
    pos = sphere1.sim.get_subfield_sites("M_Py")
    assert len(dummy) == n_nodes


    



