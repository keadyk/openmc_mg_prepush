#!/bin/bash

echo "Running SO MANY spectral radius calculations..."


#grids
for seed1 in $(seq 1 2)
do

python make_grid_file.py $seed1

#rel fact
for seed2 in $(seq 1 2)
do

python write_rel_fact.py $seed2

#sigT1
for seed3 in $(seq 1 2)
do

#sigT2
for seed4 in $(seq 1 2)
do
#sigS11
for seed5 in $(seq 1 2)
do
#sigS12
for seed6 in $(seq 1 2)
do
#sigS22
for seed7 in $(seq 1 2)
do
echo -e "RUN: $seed1 $seed2 $seed3 $seed4 $seed5 $seed6 $seed7"
python change_xs.py $seed1 $seed2 $seed3 $seed4 $seed5 $seed6 $seed7
python calc_spectral_radius_implicit_scattering_2g_relaxed_noprint.py 2g_test_input.inp 

done

done

done

done

done

done

done

echo "                              ...Finished!"
