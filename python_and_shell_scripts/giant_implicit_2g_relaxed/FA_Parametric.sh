#!/bin/bash

echo "Running SO MANY spectral radius calculations..."


#grids
for seed1 in $(seq 1 4)
do

FNAME=$(python make_grid_file.py $seed1)

#rel fact
for seed2 in $(seq 1 9)
do

#Calculate this relax factor, and write to file
RELF=$((60+5*($seed2-1))) 
#If it's 10 (maps to 1.0), we don't want the leading zero
if [ $RELF -eq 100 ] 
then
    printf "1.0 " >> "$FNAME" 
else
    printf "0.%s " "$RELF" >> "$FNAME" 
fi

#sigT1
for seed3 in $(seq 1 4)
do
#sigT2
for seed4 in $(seq 1 4)
do
#sigS11
for seed5 in $(seq 1 3)
do
#sigS12
for seed6 in $(seq 1 3)
do
#sigS22
for seed7 in $(seq 1 3)
do

#Print all zillion indices for the run
echo -e "RUN: $seed1 $seed2 $seed3 $seed4 $seed5 $seed6 $seed7"

#Update the input file
python change_xs.py $seed1 $seed2 $seed3 $seed4 $seed5 $seed6 $seed7

#Run the case, stashing the spectral radius in a shell variable
SPECRAD=$(python run_implicit_scattering_2_group_relaxed.py) 

#Append the spectral radius to the file for this grid size
#Use printf to avoid echo's newline (DUDE I DIDN'T EVEN KNOW YOU COULD DO THIS)
printf "%s " "$SPECRAD" >> "$FNAME"

#end sigS22
done

#end sigS12
done

#end sigS11
done

#end sigT2
done

#end sigT1
done

printf "\n" >> "$FNAME"
#end rel factor
done

#end grids
done

echo "                              ...Finished!"
