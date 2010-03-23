all:	run postprocess plot iterations

run:
	cp ../example2/bar30_30_100.nmesh.h5 .
	../../../../bin/nsim bar_fast.py --clean
	../../../../bin/nsim bar_normal.py --clean
	../../../../bin/nsim bar_slow.py --clean


postprocess:
	../../../../bin/ncol bar_fast time M_Py_0 M_Py_1 M_Py_2  > data_fast.txt
	../../../../bin/ncol bar_normal time M_Py_0 M_Py_1 M_Py_2  > data_normal.txt
	../../../../bin/ncol bar_slow time M_Py_0 M_Py_1 M_Py_2  > data_slow.txt

plot:
	xmgrace -nxy data_fast.txt -nxy data_normal.txt -nxy data_slow.txt -param grace.par -printfile plot.png -hardcopy
	mogrify -trim -rotate 90 plot.png

iterations:
	echo "ncol bar_fast step | tail -1" > iterations.txt
	../../../../bin/ncol bar_fast step | tail -1 >> iterations.txt
	echo "ncol bar_normal step | tail -1" >> iterations.txt
	../../../../bin/ncol bar_normal step | tail -1 >> iterations.txt
	echo "ncol bar_slow step | tail -1" >> iterations.txt
	../../../../bin/ncol bar_slow step | tail -1 >> iterations.txt
