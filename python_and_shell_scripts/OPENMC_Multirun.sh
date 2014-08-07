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
/home/keadyk/Desktop/Research_stuff/OPENMC/src/openmc &>$thisout #pipes std. out to specified output file
#move tallies.out file created by executable:
fname="${seed}_tallies.out"
mv "tallies.out" "$fname"
#echo "Moved tallies.out to $fname"

done
echo "                              ...Finished!"
#Now that you've written EVERYTHING to a CRAP TON OF FILES,
#Let the python thinger handle all the RSD calc'ing
echo "Running Python RSD calculations..."
python OPENMC_MC_RSD_calc.py $1
python OPENMC_CMFD_RSD_calc.py $1
echo "                              ...Finished!"


