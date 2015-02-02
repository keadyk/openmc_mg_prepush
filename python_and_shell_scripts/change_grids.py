import sys
import os
import math

this_run = int(sys.argv[1]) #run id-- passed from bash script

edit_file = open('test_input.inp', 'r')
temp_file = open('temp', 'w')

coarse_gridset = [1.0, 2.0, 4.0, 5.0, 6.25, 10.0]
this_gridsize = coarse_gridset[this_run-1]

for line in edit_file:
    stripped = line.strip()   #kill spaces
    if(stripped.startswith("h_j = ")):
        temp_file.write("h_j = " + str(this_gridsize) + "\n")
    else:
        temp_file.write(line)
        
        
edit_file.close()
temp_file.close()

os.rename('temp', 'test_input.inp')
