# Test the Zhang-Li extension of the LLG equation to model spin torque
import os, math
org_dir = os.getcwd()

def chdir_and_test(test):
      import sys, os

      proper_sysargv = sys.argv
      sys.argv = sys.argv[0:1]

      # work out in which directory the data files are
      org_dir = os.getcwd()
      script_dir = os.path.split(__file__)[0]
      if script_dir != "":
          os.chdir(script_dir)

      # get the data in place
      test()

      os.chdir(org_dir)


def test_run():
    import os
    def run():
        os.system("make clean m_of_t.dat > test_slow_zhangli.log")
    chdir_and_test(run)

def compare_files(file1, file2):
    """Compare two files assuming they contain a plain table of numbers
    (every line contain the same number of floating point numbers separated
    by spaces in the form Gnuplot accepts). The function checks that the two
    files contain the same number of lines and the same number of columns
    per each line. Then the corresponding entries in the two files are
    compared. The function fails when the difference is greater than the
    provided tolerance. Tolerances are specified in the list 'tolerances'.
    tolerances[i] is the tolerance for the i-th column. tolerances[-1]
    is the tolerance for the i-th column, for i >= n, where
    n = len(tolerances).
    """
    fd = open(file1, "r")
    lines1 = fd.readlines()
    fd.close()
    fd = open(file2, "r")
    lines2 = fd.readlines()
    fd.close()
    squares = []
    assert(len(lines1) == len(lines2))
    for l in range(len(lines1)):
        c1 = lines1[l].split()
        c2 = lines2[l].split()
        assert(len(c1) == len(c2))
        for i in range(len(c1)):
            n1 = float(c1[i])
            n2 = float(c2[i])
            dn = abs(n1 - n2)
	    if i >= len(squares): squares.append([0.0, 0])
	    squares[i][0] += dn*dn
	    squares[i][1] += 1
    return [math.sqrt(sum_sq_dn/num_sq_dn)
            for sum_sq_dn, num_sq_dn in squares]

def test_compare():
    def compare():
	tols = compare_files("m_of_t.dat", "m_of_t.dat.old")
        print "The differences between the two curves are", tols
	dt, dmx, dmy, dmz = tols
	assert dt < 1e-12
	assert dmx < 0.005e6
	assert dmy < 0.005e6
	assert dmz < 0.005e6
    chdir_and_test(compare)

if __name__ == "__main__":
    test_run()
    test_compare()

