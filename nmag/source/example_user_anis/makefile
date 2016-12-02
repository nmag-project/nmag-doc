.PHONY: all clean

NETGEN=netgen
GNUPLOT=gnuplot

NMESHIMPORT=../../../../bin/nmeshimport
NSIM=../../../../bin/nsim
NCOL=../../../../bin/ncol

all: coin.png

coin.neu: coin.geo
	$(NETGEN) -geofile=coin.geo -fine -meshfiletype="Neutral Format" -meshfile=coin.neu -batchmode

coin.nmesh.h5: coin.neu
	$(NMESHIMPORT) --netgen coin.neu coin.nmesh.h5

coin_dat.h5: coin.py coin.nmesh.h5
	$(NSIM) coin.py --clean

coin_hloop.dat: coin_dat.h5
	$(NCOL) coin H_ext_0 H_ext_1 H_ext_2 M_MyMat_0 M_MyMat_1 M_MyMat_2 > coin_hloop.dat

coin.eps: coin.gnp coin_hloop.dat
	$(GNUPLOT) coin.gnp

coin.png: coin.eps
	convert -trim -antialias -density 150 coin.eps coin.png
	python ../tools/scalepng coin.png 15

clean:
	rm -rf coin_dat.h5

