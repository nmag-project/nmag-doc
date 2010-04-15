NSIM_ROOT=../../src
NSIM_BIN=$(NSIM_ROOT)/bin
NSIM_PYTHONPATH=$(NSIM_ROOT)/interface
NSIM=$(NSIM_BIN)/nsim-raw

.PHONY: all pngformulae html examples latex updatefiles

all: update pngformulae html pdf updatefiles
	echo "Done."

final: update pngformulae htmlfinal pdffinal
	echo "Done."

dyn/svnversion.txt:
	echo ".. |svnversion| replace:: "$$(svnversion -n) > dyn/svnversion.txt

#the file svnversion.txt will be part of the repository for released versions
#(so we don't overwrite it with current svnversion number)
dynamic: dyn/svnversion.txt
	cd dyn; ../$(NSIM) make_dyn.py; cd ..
	cp $(NSIM_PYTHONPATH)/nsim/si_units/si.py dyn/si.py

#developers can force update of manual with proper svnversion number
#with this (non-default) target.
dev-dynamic: 
	rm dyn/svnversion.txt
	make dynamic

update: dynamic
	echo "update done."

updatefiles:
	#copy mesh files from ../examples into the right subdirectories
	cd example_ipython; make updatefiles
	cd example_larmor; make updatefiles	
	cd example_manipulate; make updatefiles
	cd example_pinning; make updatefiles
	#cd example_thermal_cube; make updatefiles

html:   
	mkdir -p html
	rst2html --strict --stylesheet-path=html/voidspace_2008_03.css --link-stylesheet -v manual.txt html/manual.html -r 1 -g -t -d -s

htmlfinal: 
	echo "Will ignore warnigs (for example 'duplicate implicite targets defined' or so)"
	mkdir -p html
	rst2html --stylesheet-path=html/voidspace.css -v manual.txt html/manual.html -r 2

latex: pdf
	echo "Done latex."

pdf: pngformulae
	rst2latex --strict --stylesheet=stylesheet.tex -v manual.txt tex/manual.tex -g -t -d
	cd tex; pdflatex manual.tex; mv -v manual.pdf ..

pdffinal: pngformulae
	echo "Will ignore warnigs (for example 'duplicate implicite targets defined' or so)"
	rst2latex --strict --stylesheet=stylesheet.tex -v manual.txt tex/manual.tex -r 2
	cd tex; pdflatex manual.tex; mv -v manual.pdf ..

htmlnmag: 
	echo "Writing dyn/nmag.txt to /tmp/test.html"
	rst2html --strict --stylesheet-path=html/voidspace.css -v dyn/nmag.txt /tmp/test.html

pngformulae:
	cd example_stoner_wohlfarth; make pngformulae
	cd example_larmor; make pngformulae
	cd tutorial; make pngformulae	

scalepngs:
	python tools/scalepng images/smallsphere.png 5

#Running examples should not be necessary as all resulting 
#image files are part of svn. It can be used as a unit test.
examples:
	cd example2; make all
	cd example3; make all
	cd example4; make all
	cd example_uniaxial_anis; make all  
	cd example_cubic_anis; make all
	cd example_pinning; make all
	cd example_thermal; make all
	cd example_1Dperiodicity; make all
	cd example_larmor; make all
	cd example_zhangli1; make all
	cd example_zhangli2; make all
