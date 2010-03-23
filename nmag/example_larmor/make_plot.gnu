set term png giant enhanced size 800,600
f(x) = A*sin(2*pi*x/B + C) + D
B = 30 
fit f(x) "data.txt" u  ($1/1e-12):2 via A,B,C,D
set out 'larmor_plot.png'
set xlabel 'Simulation time'
set ylabel 'x component of normalised magnetisation'
plot "data.txt" u ($1/1e-12):2, f(x)
