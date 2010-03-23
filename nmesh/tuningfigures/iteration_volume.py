from __future__ import division
import nmesh, pylab, timing

#create objects
box = nmesh.box([-2.0,6.0],[2.0,7.0])
ellipsoid = nmesh.ellipsoid ([1.5,1.5])

#create bounding box
bbox =([-3.,-3.],[3.,9.])

#fixed points *to fix points uncomment line below*
fix =([-2.0,6.0],[2.0,7.0],[-2.0,7.0],[2.0,6.0],
      [-3.,-3],[-3.,9.],[3.,9.],[3.,-3.])

#density
density = "density=1.0;"

#create empty string
time = []
nodes = []
average_quality=[]


#create quality stacks to plot later
stack1=[]
stack2=[]
stack3=[]
stack4=[]
stack5=[]
stack6=[]

iteration = [1,5,10,50,200,500,1000,2000]

#mesh object at different maximum number of iterations N
for N in iteration:
    timing.start()
    mesh = nmesh.mesh(objects=[box, ellipsoid], bounding_box=bbox, \
                      a0=0.3, mesh_bounding_box=True,
                      
                      fixed_points = fix,
 
                      max_steps = N,

                      density = density,

                      neigh_force_scale = 0.0,

                      shape_force_scale = 0.0,

                      volume_force_scale = 1.0
                      )
    timing.finish()
    time.append(timing.seconds())
    #iteration.append(N)

    #save the mesh as a .ps file in temp dir
    nmesh.visual.plot2d_ps(mesh,"fig_mesh_volume_iter%06d.ps" % N) 

    #extract the mesh data from mesh_info
    mesh_info = mesh.tolists()
    vtkData, points, simplices, simplexIndicies, icradii, ccradii = \
             nmesh.visual.mesh2vtk(mesh_info, VTKonly=False)
    in2circ = nmesh.visual.findRatios(icradii, ccradii, factor=2)#2D

    #count the number of point after iterations and add to list
    nodes.append(len(points))

    #calculate the average ratio and append to list
    quality= sum(in2circ)/len(in2circ)
    average_quality.append(quality)

    #create counters
    counter1 = 0
    counter2 = 0
    counter3 = 0
    counter4 = 0
    counter5 = 0
    counter6 = 0
    
    #iterate through the list of ratio and obtain the intervals counts
    for i in range(len(in2circ)):
        if in2circ[i]<=0.1:
            counter1 = counter1 + 1
        else:
            counter1 = counter1

        
        if in2circ[i]<=0.25 and in2circ[i]>0.1:
            counter2 = counter2 + 1
        else:
            counter2 = counter2
        

        if in2circ[i]<=0.50 and in2circ[i]>0.25:
            counter3 = counter3 + 1
        else:
            counter3 = counter3
        
    
        if in2circ[i]<=0.75 and in2circ[i]>0.50:
            counter4 = counter4 + 1
        else:
            counter4 = counter4
        
                
        if in2circ[i]<=0.90 and in2circ[i]>0.75:
            counter5 = counter5 + 1
        else:
            counter5= counter5
        

        if in2circ[i]<=1.0 and in2circ[i]>0.90:
            counter6 = counter6 + 1
        else:
            counter6= counter6
        

    #calculate the percentage of the simplices within each interval        
    counter1 = counter1*100/len(in2circ)
    counter2 = counter2*100/len(in2circ)+ counter1
    counter3 = counter3*100/len(in2circ)+ counter2
    counter4 = counter4*100/len(in2circ)+ counter3
    counter5 = counter5*100/len(in2circ)+ counter4
    counter6 = counter6*100/len(in2circ)+ counter5

    #append the counters to the list
    stack1.append(counter1)
    stack2.append(counter2)
    stack3.append(counter3)
    stack4.append(counter4)
    stack5.append(counter5)
    stack6.append(counter6)
    
    #plot histogram of mesh quality
    pylab.ion()
    pylab.clf
    pylab.hist(in2circ, bins=50)
    pylab.ylim(ymax=1600)
    pylab.title("Quality of mesh after %d iterations \n(Time \
    taken to mesh object = %d seconds)" % (N, timing.seconds()))
    pylab.xlabel("2*inradius/circumradius")
    pylab.ylabel("Number of simplices")
    #pylab.draw()
    pylab.hold()
    pylab.savefig("fig_hist_volume_iter%06d.eps" % N) #save histogram
    pylab.ioff
    pylab.hold(False)

#plot line graph of iteration vs time
pylab.clf
pylab.plot(iteration, time, 'o-')
pylab.title(" Time for iterations")
pylab.ylabel("Time (seconds)")
pylab.xlabel("Number of iterations")
pylab.savefig("graph_time_volume.eps")

#plot line graph of iteration vs number of points
pylab.clf
pylab.plot(iteration, nodes, 'o-')
pylab.title(" Number of points after iterations")
pylab.ylabel("Number of nodes")
pylab.xlabel("Number of iterations")
pylab.savefig("graph_node_volume.eps")

#plot line graph of iteration vs average_quality
pylab.clf
pylab.plot(iteration,average_quality)
pylab.title("Average quality")
pylab.ylabel("Average 2*inradius/circumradius ratio")
pylab.xlabel("iterations")
pylab.savefig("graph_averagequality_volume.eps")

#plot line graph of iteration vs the percentage of simplice
#in each interval
pylab.clf
pylab.plot(iteration, stack1,label='0-10%')
pylab.hold()
pylab.plot(iteration, stack2, label='10-20%')
pylab.plot(iteration, stack3, label='25-50%')
pylab.plot(iteration, stack4, label='50-75%')
pylab.plot(iteration, stack5, label='75-90%')
pylab.plot(iteration, stack6,label='90-100%')
pylab.legend()
pylab.title("Quality Stacks")
pylab.ylabel("Percentage of simplices with ratio within each interval")
pylab.xlabel("iterations")
pylab.savefig("graph_quality_volume.eps")
