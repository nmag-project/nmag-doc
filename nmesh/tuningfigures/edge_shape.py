from __future__ import division
import nmesh, pylab, timing, Numeric


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

force=[0,1,2,5,7,10,20,50,70,100] 
for e in force:
    shape_force = e/100.0
    neigh_force = 1. - shape_force

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

    iteration = [1000]


    #mesh object at different maximum number of iterations N
    for N in iteration:
        timing.start()
        mesh = nmesh.mesh(objects=[box, ellipsoid], bounding_box=bbox,\
                          a0=0.3, mesh_bounding_box=True,

                          fixed_points = fix,

                         
                          #tolerated_rel_move = 0.0009,

                          max_steps = N,

                          density = density,

                          neigh_force_scale = neigh_force,

                          shape_force_scale = shape_force,

                          volume_force_scale = 0.0
                          )
        timing.finish()
        time.append(timing.seconds())

        # save mesh
        mesh.save("edge_shape_%d.nmesh" % e)

        # load mesh from file
        #mesh_loaded = nmesh.load("save_load.nmesh")

        # plot loaded mesh
        #nmesh.visual.plot2d_ps(mesh_loaded,"save_load.ps")


        #save the mesh as a .ps file in temp dir
        nmesh.visual.plot2d_ps\
             (mesh,"fig_mesh_edge_shape_%d_iter%06d.ps" % (e,N)) 

        #extract the mesh data from mesh_info
        mesh_info = mesh.tolists()
        vtkData, points, simplices, simplexIndicies, icradii, ccradii=\
nmesh.visual.mesh2vtk(mesh_info, VTKonly=False, body_numbers=False)
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

        #compute the interval counts
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


        #calculate the percentage of the simplices in each interval        
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
        pylab.ylim(ymax=1800)
        pylab.title("Quality of mesh after %d iterations\n\
        (Time taken to mesh object = %d seconds)"\
                    % (N, timing.seconds()))
        pylab.xlabel("2*inradius/circumradius")
        pylab.ylabel("Number of simplices")
        #pylab.draw()
        pylab.hold()
        pylab.savefig("fig_hist_edge_shape_%d_iter%06d.eps"\
                      % (e, N)) #save histogram
        pylab.ioff
        pylab.hold(False)

        #get the bad simplices
        ratiobad = []
        for  i in range(len(in2circ)):
            if in2circ[i]< 0.80:
                ratiobad.append(in2circ[i])

        if len(ratiobad) != 0:
            #plot a histogram of the bad simplices
            pylab.ion()
            pylab.clf
            pylab.grid()
            pylab.hist(ratiobad, bins=50)
            pylab.ylim(ymax=20)
            pylab.xlim(xmin=0.0, xmax=0.8)
            pylab.title\
("Number of low quality  simplices after %d iterations\n\
(Time taken to mesh object = %d seconds)" % (N, timing.seconds()))
            pylab.xlabel("2*inradius/circumradius")
            pylab.ylabel("Number of simplices")
            #pylab.draw()
            pylab.hold()
            pylab.savefig\
("fig_hist_edge_shape_%d_badsimplices_iter%06d.eps"\
 %  (e, N)) #save histogram
            pylab.ioff
            pylab.hold(False)

       #similary get the good simplices 
        ratiogood = []
        for  i in range(len(in2circ)):
            if in2circ[i]>=0.80:
                ratiogood.append(in2circ[i])

        #plot a histogram of the good simplices
        pylab.ion()
        pylab.clf
        pylab.hist(ratiogood, bins=50)
        pylab.ylim(ymax=1600)
        pylab.title\
("Number of high quality simplices after %d iterations\n\
(Time taken to mesh object = %d seconds)" % (N, timing.seconds()))
        pylab.xlabel("2*inradius/circumradius")
        pylab.ylabel("Number of simplices")
        pylab.hold()
        pylab.savefig\
("fig_hist_edge_shape_%d_goodsimplices_iter%06d.eps" %  (e, N)) 
        pylab.ioff
        pylab.hold(False)

         #similary get the good simplices 
        ratiogood = []
        for  i in range(len(in2circ)):
            if in2circ[i]>=0.80:
                ratiogood.append(in2circ[i])

        #plot a histogram of the good simplices
        pylab.ion()
        pylab.clf
        pylab.hist(ratiogood, bins=50)
        pylab.ylim(ymax=1600)
        pylab.title\
("Number of high quality simplices after %d iterations\n\
(Time taken to mesh object = %d seconds)" % (N, timing.seconds()))
        pylab.xlabel("2*inradius/circumradius")
        pylab.ylabel("Number of simplices")
        pylab.hold()
        pylab.savefig\
("fig_hist_edge_shape_%d_goodsimplices_iter%06d.eps" %  (e, N)) 
        pylab.ioff
        pylab.hold(False)

