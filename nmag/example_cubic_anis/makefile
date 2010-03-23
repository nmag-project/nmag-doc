run: simulate media

simulate:
	../../../../bin/nsim cube.py --clean


ncol:
	../../../../bin/ncol cube H_ext_0 M_Fe_0 > cube_hext_vs_m.txt

media:
	gnuplot cube_hext_vs_m.gnu
	python ../tools/scalepng cube_hext_vs_m.png 15
