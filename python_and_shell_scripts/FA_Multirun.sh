#!/bin/bash

echo "Running k-eigenvalue calculations..."

for seed1 in $(seq 1 6)
do
python change_grids.py $seed1

for seed2 in $(seq 1 12)
do

echo -e "RUN: $seed1 $seed2"
python change_scatrat.py $seed2
python calc_spectral_radius_n_inner_scattering_Larsen_part_sum.py test_input.inp 

done

done

echo "                              ...Finished!"
