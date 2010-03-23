for x in *py; do 
    echo "working on $x" 
    ../../../../bin/nmesh2 $x > ${x%py}log; 
done


for x in *.ps; do
    echo "working on $x" 
    ps2epsi $x ${x%ps}eps
    rm $x
done

rm *vtk



