import pylab

#open  the file iteration.log
f=open('iteration.log','r')

#create call_counter list
iteration=[]
counter_calls=[]

#for N in [1,5,10]:#,20,30,40,50,200,500,1000,2000]

#read file line by line
readlines = f.readlines()

#set counter to zero

for N in [1,5,10,20,30,40,50,200,500,1000,2000]:
    counter = 0
    for line in readlines:
        
        if line[44:53]== "MAX_STEPS":
            bits = line.split()
            #print bits
            max_step = int( bits[5] )#convert bits[2] to an interger
            #print max_step

            if N == max_step:    #if max_step = N increase counter by 1
                counter=counter+1
                #print counter
            #else:
                #counter = counter #else counter does not change
                #print counter
    counter_calls.append(counter) #add the final value to list
    iteration.append(N)
#close the file
f.close()

print counter_calls
pylab.plot(iteration,counter_calls,'o-')
pylab.title("Calls to Qhull")
pylab.ylabel("Number of calls to Qhull")
pylab.xlabel("Iterations")
pylab.savefig("graph_qhull_call_edge.eps")

