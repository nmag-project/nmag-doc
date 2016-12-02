#!/bin/bash

if [ -e timings_full.dat ]
then
cp timings_full.dat timings_full.dat_old
fi
../../../../../bin/nsim thinfilm20_full.py --clean
../../../../../bin/nsim thinfilm40_full.py --clean
../../../../../bin/nsim thinfilm60_full.py --clean
../../../../../bin/nsim thinfilm70_full.py --clean
../../../../../bin/nsim thinfilm80_full.py --clean
../../../../../bin/nsim thinfilm90_full.py --clean
../../../../../bin/nsim thinfilm100_full.py --clean
../../../../../bin/nsim thinfilm110_full.py --clean
../../../../../bin/nsim thinfilm120_full.py --clean
../../../../../bin/nsim thinfilm130_full.py --clean
