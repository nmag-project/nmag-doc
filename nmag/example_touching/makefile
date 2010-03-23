BIN=../../../../bin
NMESHIMPORT=$(BIN)/nmeshimport
NETGEN=netgen

touching.nmesh.h5: touching.geo
	netgen -geofile=$< -meshfiletype="Neutral Format" -meshfile=mesh.neu -batchmode
	$(NMESHIMPORT) --netgen mesh.neu $@
	rm -f mesh.neu

touching.geo: touching.geo.py
	python $< $@

