#!/usr/bin/python

import sys, os

otheropts = ""
if len(sys.argv) == 3:
  _, filename, formula = sys.argv
else:
  try:
    _, filename, formula, border = sys.argv
    otheropts = "-border %s -bordercolor white" % border
  except:
    print "USAGE:  python latex-png.py outfile.png \"my = eq\" [border]"
    print "EXAMPLE 1: python latex-png.py sqrt.png \"c = \sqrt{a^2 + b^2}\""
    print "EXAMPLE 2: python latex-png.py my_equation.png \"c = \sqrt{a^2 + b^2}\" 5"
    sys.exit(1)

tmpdir = "/tmp/latex-formula-%d" % os.getpid()
os.mkdir(tmpdir)

tex = """
\\documentclass{article}
\\pagestyle{empty}
\\begin{document}
\\begin{eqnarray}
%s
\\nonumber
\\end{eqnarray}
\\end{document}
""" % formula

f = open("%s/formula.tex" % tmpdir, "wt")
f.write(tex)
f.close()

os.system("cd %s; latex formula.tex; dvips -f <formula.dvi >formula.ps; ps2epsi formula.ps formula.eps; convert -density 131 %s formula.eps formula.png" % (tmpdir, otheropts));

os.system("cp %s/formula.png %s" % (tmpdir, filename));

os.system("rm -rf %s" % tmpdir);
