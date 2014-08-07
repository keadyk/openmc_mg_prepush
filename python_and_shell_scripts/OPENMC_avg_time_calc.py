from __future__ import division
import sys
import os.path
import math

nfiles = int(sys.argv[1]) #Number of files to look for-- passed from bash script

file_string1 = "avg_time.out"

dest_file1 = open(file_string1, 'w')
print "Created output file " + file_string1
dest_file1.write("Avg. simulation time: ")

realfiles = nfiles
index = 0

time_sum = 0

for l in range(nfiles): #For the expected number of files...
    index = l + 1 #1-indexed file number
    if os.path.isfile(`index` + "_timing.out") is True:
        print "Found file " + `index` + "_timing.out!"
        this_file = open(`index` + "_timing.out", 'r')
    else:
        print "File " + `index` + "_timing.out is missing!"
        realfiles -= 1 #decrement the number of 'real' files
        continue
    for line in this_file:
        #split into values by whitespace
        values = line.split()
        if(len(values) == 0): 
            continue
        elif(values[0].startswith("Reading")):
            continue
        #If we make it to here, it isn't a blank or comment line!
        elif(values[0].startswith("Total")):
            if(values[1].startswith("time")):
                if(values[2].startswith("elapsed")):
                    #We found the total time elapsed! Add to sum...
                    time_sum += float(values[4])

#All the error lines should be the same length
time_avg = time_sum/float(realfiles)

dest_file1.write(str(time_avg))

dest_file1.close()

print "Done with avg time calc.  See " + file_string1 + "."  




        
