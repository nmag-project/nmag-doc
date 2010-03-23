#!/bin/bash

if [ -e timings_hlib.dat ]
then
cp timings_hlib.dat timings_hlib.dat_old
fi
../../../../../bin/nsim thinfilm20_hlib.py --clean
../../../../../bin/nsim thinfilm40_hlib.py --clean
../../../../../bin/nsim thinfilm60_hlib.py --clean
../../../../../bin/nsim thinfilm70_hlib.py --clean
../../../../../bin/nsim thinfilm80_hlib.py --clean
../../../../../bin/nsim thinfilm90_hlib.py --clean
../../../../../bin/nsim thinfilm100_hlib.py --clean
../../../../../bin/nsim thinfilm110_hlib.py --clean
../../../../../bin/nsim thinfilm120_hlib.py --clean
../../../../../bin/nsim thinfilm130_hlib.py --clean
