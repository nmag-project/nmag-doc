geo = """
algebraic3d

solid cube1 = orthobrick (    0, -$L$, -$L$; $L2$, $L$, $L$) -maxh = $D$;
solid cube2 = orthobrick (-$L2$, -$L$, -$L$;    0, $L$, $L$) -maxh = $D$;
tlo cube1;
tlo cube2;
"""

L = 10
D = 5
var_vals = [("L", L), ("L2", 2*L), ("D", D)]

import sys
f = open(sys.argv[1], "w")
for var, val in var_vals: geo = geo.replace("$%s$" % var, str(val))
f.write(geo)
f.close()

