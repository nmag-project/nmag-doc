# I'm using Python merely to detect the machine on which I am running,
# and to insert the relevant flags for that platform/directory
# structure. It's probably easy to do within the Makefile itself
# (which would be preferable), but I'd rather get on with the problem
# at hand, rather than learning about what is possible in Makefiles.
######################################################################

######################################################################
eta_opts = \
"""
PY_PREFIX=/home/jacek
PY_VERSION=2.5

LIBS = -L/home/jacek/three_eleven/nmag-0.1/lib/ocaml/ -L$(PY_PREFIX)/lib/python$(PY_VERSION)/config -lpython$(PY_VERSION) -lunix -lm -lcurses 
INCS = -I$(PY_PREFIX)/include/python$(PY_VERSION) -I/home/jacek/three_eleven/nmag-0.1/lib/ocaml/

CC_SHARED_FLAGS = -shared --whole-archive 
"""
######################################################################
orwell_opts = \
"""
PY_PREFIX=/usr/local
PY_VERSION=2.4

LIBS = -L/usr/local/lib/ocaml/ -L$(PY_PREFIX)/lib/python$(PY_VERSION)/config -lpython$(PY_VERSION) -lunix -lm -lcurses 
INCS = -I$(PY_PREFIX)/include/python$(PY_VERSION) -I/usr/local/lib/ocaml/

CC_SHARED_FLAGS = -bundle -flat_namespace -undefined suppress
"""
######################################################################
targets = \
"""
all:
	make clean
	make byte
	make testbyte
	make testbyte2
	make native
	make test
	make test2


native: themodule.so

themodule.o: themodule.c
	gcc -ggdb -Wall -c themodule.c $(INCS) -o themodule.o

mlcore.o: mlcore.ml
	ocamlopt -o mlcore.so -output-obj -ccopt -fPIC mlcore.ml

themodule.so: themodule.o mlcore.o
	gcc $(CC_SHARED_FLAGS) mlcore.so themodule.o $(LIBS) -lasmrun -o themodule.so

test: themodule.so
	python$(PY_VERSION) test.py

test2: themodule.so
	python$(PY_VERSION) -c "import themodule; print 'square(2) =', themodule.pysquare(2)"

testnative: test

testnative2: test2




byte: thebytemodule.so

thebytemodule.o: thebytemodule.c
	gcc -ggdb -Wall -c thebytemodule.c $(INCS) -o thebytemodule.o

bytemlcore.o: mlcore.ml
	ocamlc -o bytemlcore.o -output-obj -ccopt -fPIC mlcore.ml

thebytemodule.so: thebytemodule.o  bytemlcore.o
	gcc $(CC_SHARED_FLAGS)  bytemlcore.o thebytemodule.o $(LIBS) -lcamlrun -o thebytemodule.so

testbyte: thebytemodule.so
	python$(PY_VERSION) testbyte.py

testbyte2: thebytemodule.so
	python$(PY_VERSION) -c "import thebytemodule; print 'square(2) =', thebytemodule.pysquare(2)"

clean:
	rm -rf *o *so *cmi *cmo *cmx
"""
######################################################################
import os
if os.uname()[0] == 'Linux':
    opts = eta_opts
else:
    opts =  orwell_opts
makefile = open('Makefile','w')
print >> makefile, opts
print >> makefile, targets
