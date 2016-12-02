lines = open('resultsummary.txt','r').readlines()

data_by_tol = {}

for line in lines:
    data = line.split()
    data_by_tol[float(data[0])] = int(data[1]),float(data[2])

tols = data_by_tol.keys()

tols.sort()

f=open('resultsummary_rst.txt','w')

f.write('%10s %10s %14s %21s\n' % ('==========','==========','==============','====================='))
f.write('%10s %10s %14s %21s\n' % ('tol','steps','CPU time (s)','CPU time per step (s)'))
f.write('%10s %10s %14s %21s\n' % ('==========','==========','==============','====================='))


for tol in tols:
    steps,cpu_time = data_by_tol[tol]
    cpu_time_per_step = cpu_time / steps
    f.write('%10.6f %10d %14.2f %21.3f\n' % ( tol, steps,cpu_time, cpu_time_per_step ))

f.write('%10s %10s %14s %21s\n' % ('==========','==========','==============','====================='))

f.close()
