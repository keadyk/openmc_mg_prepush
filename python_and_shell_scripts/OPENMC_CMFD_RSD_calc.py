from __future__ import division
import sys
import os.path
import math

nfiles = int(sys.argv[1]) #Number of files to look for-- passed from bash script

file_string1 = "CMFD_Real_RSD.out"

dest_file1 = open(file_string1, 'w')
print "Created output file " + file_string1
dest_file1.write("Mesh indices, Avg. Flux, RSD: ")

realfiles = nfiles
index = 0
fileflag = 1
Cindices = [[],[],[]]

Cval = []
Cvalsq = []

tallyid = 0

for l in range(nfiles): #For the expected number of files...
    index = l + 1 #1-indexed file number
    if os.path.isfile(`index` + "_CMFD.out") is True:
        print "Found file " + `index` + "_CMFD.out!"
        this_file = open(`index` + "_CMFD.out", 'r')
        i = 0
    elif os.path.isfile(`index` + "_HCMFD.out") is True:
        print "Found file " + `index` + "_HCMFD.out!"
        this_file = open(`index` + "_HCMFD.out", 'r')
        i = 0
    else:
        print "Index # " + `index` + " is missing!"
        realfiles -= 1 #decrement the number of 'real' files
        continue
    for line in this_file:
        #split into values by whitespace
        values = line.split()
        if(len(values) == 0): 
            continue
        elif(values[0].startswith("#")):
            continue
        #If we make it to here, it isn't a blank or comment line!
        if(fileflag == 1): #This is the first file we've looked at
            #print "Getting indices " + values[2].strip('(,)') + " " + values[3].strip('(,)') + " " + values[4].strip('(,)')
            #strip all commas, parens. Grab indices!
            j = int(values[0])
            k = int(values[1])
            m = int(values[2])
            #store indices:
            Cindices[0].append(j)
            Cindices[1].append(k)
            Cindices[2].append(m)			        
            #store values!
            Cval.append(float(values[4]))
            Cvalsq.append(float(values[4]) * float(values[4]))
        else:
            Cval[i] += (float(values[4]))
            Cvalsq[i] += (float(values[4]) * float(values[4]))
        i += 1 #Increment the line counter
        #print i
    #unset the first-file flag
    fileflag = 0

Ccells = len(Cindices[0])

CRSD = [0]*Ccells
	
for i in range(Ccells):
	#Calc the RSD value 
	CRSD[i] = math.fabs( math.sqrt((realfiles/(realfiles - 1.0)) * math.fabs((Cvalsq[i])/realfiles - ((Cval[i])*(Cval[i]))/(realfiles*realfiles) ) ) / (Cval[i]/realfiles) )
	#print "RSD [" + str(g) + "][" + str(i) + "]: " + str(RSD[g][i])

#All the error lines should be the same length
for i in range(Ccells):
	this_string = str(Cindices[0][i]) + '\t' + str(Cindices[1][i]) + '\t' + str(Cindices[2][i]) + '\t' + str(Cval[i]/realfiles) + '\t' + str(CRSD[i]) #write cell indices, RSD
	dest_file1.write('\n' + this_string)

dest_file1.close()

print "Done with real RSD calc.  See " + file_string1 + " and for real RSD data."  




        
