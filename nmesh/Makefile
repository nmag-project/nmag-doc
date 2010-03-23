all:
	make figure
	make tex


figure:
	cd figures; sh update_eps.sh; cd ..

tex:
	latex nmesh.tex
	latex nmesh.tex
	dvips -Ppdf -onmesh.ps nmesh.dvi
	ps2pdf nmesh.ps



#Kon's work (Sep 2006)
#tuningfigure:
#	cd tuningfigures; sh update_eps.sh; cd ..
#
#
#tuningtex:
#	latex tuning.tex
#	latex tuning.tex
#	dvips -Ppdf -otuning.ps tuning.dvi
#	ps2pdf tuning.ps
#
