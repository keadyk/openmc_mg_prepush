#!/bin/bash

if [ $# != 1 ]
then
echo "Usage: bash OPENMC_MultiMultirun.sh <# RUNS>"
exit 1
fi

echo "Number of runs for each simulation: $1"
echo "Running FIRST real RSD simulation..."

cd ./34_34_x5x5focus_real_rsd
/bin/bash ./OPENMC_Multirun.sh $1

echo "Finished first real RSD simulation..."
echo "Starting second real RSD simulation..."

cd ../34_34_x25focus_real_rsd
/bin/bash ./OPENMC_Multirun.sh $1
