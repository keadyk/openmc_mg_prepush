import sys
import os
import math

this_run = int(sys.argv[1]) #run id-- passed from bash script

edit_file = open('test_input.inp', 'r')
temp_file = open('temp', 'w')

srat_set = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.99]
this_srat = srat_set[this_run-1]

for line in edit_file:
    stripped = line.strip()   #kill spaces
    if(stripped.startswith("s_rat = ")):
        temp_file.write("s_rat = " + str(this_srat) + "\n")
    else:
        temp_file.write(line)
        
        
edit_file.close()
temp_file.close()

os.rename('temp', 'test_input.inp')
