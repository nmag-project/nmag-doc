"""This test program is called by py.test (from many directory levels above)."""

import os

def test_slow_hysteresis_ellipsoid():
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

    #compute new data (this can take a while, will delete old data files automatically)
    os.system("make run > test_make_run1.out")

    import nmag

    #load saved magnetisation before switching

    before_switching=nmag.get_subfield_from_h5file('ellipsoid_dat.h5','m_Py',id=26)

    Mx = before_switching[:,0]
    My = before_switching[:,1]
    Mz = before_switching[:,2]

    average_Mx = Mx.sum()/len(Mx)
    average_My = My.sum()/len(Mx)
    average_Mz = Mz.sum()/len(Mx)

    #with svnrev=5886: Average mx, my, mz = (0.995631 -0.093231 0.002322)
    print "Average mx, my, mz = (%f %f %f)" % (average_Mx,average_My,average_Mz)

    assert average_Mx > 0.99, "Expect <Mx> at id=26 to be > 0.99 but is %f" % average_Mx
    assert abs(average_My) < 0.1, "Expect <My> at id=26 to be < 0.1 but is %f" % average_My
    assert abs(average_Mz) < 0.01, "Expect <Mz> at id=26 to be < 0.1 but is %f" % average_Mz


    after_switching=nmag.get_subfield_from_h5file('ellipsoid_dat.h5','m_Py',id=27)

    Mx = after_switching[:,0]
    My = after_switching[:,1]
    Mz = after_switching[:,2]

    average_Mx = Mx.sum()/len(Mx)
    average_My = My.sum()/len(Mx)
    average_Mz = Mz.sum()/len(Mx)

    #with svnrev=5886: Average mx, my, mz = (-0.999987 -0.005147 -0.000053) 
    print "Average mx, my, mz = (%f %f %f)" % (average_Mx,average_My,average_Mz)

    assert average_Mx < -0.99, "Expect <Mx> at id=26 to be > 0.99 but is %f" % average_Mx
    assert abs(average_My) < 0.1, "Expect <My> at id=26 to be < 0.1 but is %f" % average_My
    assert abs(average_Mz) < 0.01, "Expect <Mz> at id=26 to be < 0.1 but is %f" % average_Mz

    #this avoids getting relative paths in the checking progress output
    os.chdir(org_dir)



