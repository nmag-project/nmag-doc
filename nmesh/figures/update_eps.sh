FILELISTRUN="tutorial1 tutorial2 tutorial3 box frustum frustum3d shift scale rotate2d rotate3d rotate transformations multiobjects union difference intersection save_load simple1d simple2d simple3d simple4d ellipsoid_array ellipsoid_array3d xellipsoid_array3d mergedspheres periodic density1 density2 fixedpoints helix xellipsoid_array3d"

FILELISTFIGURES="tutorial1 tutorial2 tutorial3 box frustum frustum3d shift scale rotate2d rotate2d rotate transformations multiobjects union difference intersection save_load simple2d simple3d ellipsoid_array ellipsoid_array3d  mergedspheres periodic density1 density2 fixedpoints"


for x in $FILELISTRUN; do 
    echo "working on $x" 
    echo "working on $x" >> /tmp/progress.txt
    ../../../../bin/nsim ../../examples/${x}.py  ; 
done

for x in $FILELISTFIGURES; do
    pushd .
    cd run_${x%.py}
    echo "working on ${x%.py}" 
    #ls -l
    for y in *.ps; do
	ps2epsi ${y} ${y%.ps}.eps
    done
    #copy into figures directory
    cp *.eps ../
    popd
done


cp run_rotate3d/rotate3d.eps .

cp run_helix/helix.eps .

cp run_simple1d/simple1d.nmesh .

cp run_simple3d/simple3d.eps .

cp run_ellipsoid_array3d/ellipsoid_array3d.eps .

gunzip -f *eps.gz


