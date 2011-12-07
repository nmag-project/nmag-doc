""".. |nmagverstr| replace:: 0.1.0-dev
.. |nmagverinfo| replace:: 0.1.0-dev (?) released on ?"""
import sys
from nsim.version import version_str, release_date, vcinfo

if "--info" in sys.argv:
    rd = (" released on %s" % release_date) if release_date else ""
    vc = (" (%s)" % (vcinfo.split()[1])) if vcinfo else ""
    print version_str + vc + rd

else:
    print version_str

