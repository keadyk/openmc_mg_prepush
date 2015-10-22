#!/bin/bash

if [ $# != 1 ]
then
echo "Usage: bash OPENMC_Multirun.sh <# RUNS>"
exit 1
fi

echo "Number of runs: $1"
echo "Running k-eigenvalue calculations..."

for seed in $(seq 1 $1)
do

echo -e "RUN: $seed"
thisout="${seed}_timing.out"
python OPENMC_chg_seed.py $seed
/home/keadyk/OpenMC/src/openmc_SPECRAD &>$thisout #pipes std. out to specified output file
#move tallies.out file created by executable:
fname="${seed}_tallies.out"
mv "tallies.out" "$fname"
#echo "Moved tallies.out to $fname"

done
echo "                              ...Finished!"

