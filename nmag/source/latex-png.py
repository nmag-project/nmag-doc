#!/usr/bin/python

def usage():
    print "USAGE:  python latex-png.py outfile.png \"my = eq\" [border]"
    print "EXAMPLE 1: python latex-png.py sqrt.png \"c = \sqrt{a^2 + b^2}\""
    print "EXAMPLE 2: python latex-png.py my_equation.png \"c = \sqrt{a^2 + b^2}\" 5"
    print "EXAMPLE 2: python latex-png.py -f equation.tex my_equation.png 5"
    sys.exit(1)

def get(l, i, default=None):
  if i < len(l):
    return l[i]
  elif default == None:
    usage()
  else:
    return default

import sys, os

args = sys.argv[1:]
l = len(args)
if l < 1:
  usage()

border = False
if args[0] == "-f":
  assert l >= 2, \
    "-f must be followed by the name of the file containing the equation"
  f = open(args[1], "r")
  formula = f.read().strip()
  f.close()

  filename = get(args, 2)
  border = get(args, 3, False)

elif l == 2:
  filename, formula = args

elif l == 3:
  filename, formula, border = args

otheropts = "-border %s -bordercolor white" % border if border != False else ""
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
