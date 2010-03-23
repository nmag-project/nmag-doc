"""This test program is called by py.test (from many directory levels above)."""

import os

#delete old data files to avoid failure due to this
def remove_if_exists(filename):
    if os.path.exists(filename):
        os.remove(filename)


def test_slow_mfm():
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

    #remove old data file
    remove_if_exists('disk_dat.h5')

    #get the mesh in place (this requires netgen, so we have
    #included the mesh file in the repository)
    #os.system("make remesh > test_make_remesh.out")

    #compute new data (this can take a while)
    os.system("make run1 > test_make_run1.out")

    import nmag

    #load saved magnetisation

    relaxed_m=nmag.get_subfield_from_h5file('disk_dat.h5','m_Py',row=0)

    Mx = relaxed_m[:,0]
    My = relaxed_m[:,1]
    Mz = relaxed_m[:,2]

    average_Mx = Mx.sum()/len(Mx)
    average_My = My.sum()/len(Mx)
    average_Mz = Mz.sum()/len(Mx)

    print "Average mx, my, mz = (%f %f %f)" % (average_Mx,average_My,average_Mz)

    #we expect a vortex, so Mx and My should be fairly small, and Mz non-zero
    #On 23 July 2008 (rev 5640), I get
    #   """Average mx, my, mz = (0.001327 0.006164 0.031480)"""
    #Let's check for that:

    assert abs(average_Mx) < 0.01, "Mx=%f is too large" % average_Mx
    assert abs(average_My) < 0.01, "My=%f is too large" % average_My
    assert abs(average_Mz) > 0.03, "Mz=%f is too small" % average_Mz

    #this avoids getting relative paths in the checking progress output
    os.chdir(org_dir)


