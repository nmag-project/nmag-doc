import sys
from nsim.version import version_str, release_date, vcinfo

with open(sys.argv[1], "r") as f:
    content = f.read()

rd = (" released on %s" % release_date) if release_date else ""
vc = (" (%s)" % (vcinfo.split()[1])) if vcinfo else ""

verstr = (version_str[:-4] if version_str.endswith("-dev") else version_str)
verinfo = verstr + vc + rd
content = (content.replace("@NMAG_VERSION@", verstr)
                  .replace("@NMAG_VERINFO@", verinfo))

with open(sys.argv[2], "w") as f:
    f.write(content)
