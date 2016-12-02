SRC=../example_cubic_anis/cube.py
SRC_MESH=../example_cubic_anis/cube.nmesh
DST=cube_restartable.py
NSIM=../../../../bin/nsim
NCOL=../../../../bin/ncol

.PHONY: all run restart

all: $(DST) last_line.py

$(DST): $(SRC)
	NL=`cat $(SRC) | wc -l`; \
	 NLm1=`echo $$NL-1 | bc`; \
	 head -n$$NLm1 $(SRC) > $(DST) && \
	 cat ./restart_modifications.py >> $(DST)

last_line.py: $(SRC)
	tail -n1 $(SRC) > last_line.py

cube.nmesh:
	ln -s $(SRC_MESH) .

run: cube.nmesh $(DST)
	$(NSIM) $(DST) --clean

restart: cube.nmesh $(DST)
	$(NSIM) $(DST) --restart
