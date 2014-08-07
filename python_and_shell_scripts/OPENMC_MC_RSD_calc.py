from __future__ import division
import sys
import os.path
import math

nfiles = int(sys.argv[1]) #Number of files to look for-- passed from bash script

file_string1 = "MC_Real_RSD.out"
file_string2 = "Tally2_Real_RSD.out"

dest_file1 = open(file_string1, 'w')
dest_file2 = open(file_string2, 'w')
print "Created output files " + file_string1 + " and " + file_string2
dest_file1.write("Mesh indices, RSD: ")
dest_file2.write("Mesh indices, RSD: ")

realfiles = nfiles
index = 0
fileflag = 1
Findices = [[],[],[]]
Cindices = [[],[],[]]

Fval = []
Fvalsq = []

Cval = []
Cvalsq = []

tallyid = 0

for l in range(nfiles): #For the expected number of files...
    index = l + 1 #1-indexed file number
    if os.path.isfile(`index` + "_tallies.out") is True:
        print "Found file " + `index` + "_tallies.out!"
        this_file = open(`index` + "_tallies.out", 'r')
        i = 0
        for line in this_file:
            #split into values by whitespace
            values = line.split()
            if(len(values) == 0): 
                continue
            elif(values[0].startswith("===")):
                #This is a tally identifier line!
                if(int(values[2].strip(':')) == 1):
                    tallyid = 1
                    continue  
                elif(int(values[2].strip(':')) == 2):
                    #!!!IMPORTANT!!! Reset line counter:
                    i = 0
                    tallyid = 2
                    continue
                else:
                    break   #This is not a tally we want!
            elif(values[0].startswith("Total") or values[0].startswith("Vol")):
                continue    #skip this line
            if(fileflag == 1): #First file we've looked at!
                if(values[0].startswith("Mesh")):
                    #print "Getting indices " + values[2].strip('(,)') + " " + values[3].strip('(,)') + " " + values[4].strip('(,)')
                    #strip all commas, parens. Grab indices!
                    j = int(values[2].strip('(,)'))
                    k = int(values[3].strip('(,)'))
                    if(len(values) > 4):
                        m = int(values[4].strip('(,)'))
                    else:
                        m = 1
                    #store indices:
                    if(tallyid == 1):
                        Findices[0].append(j)
                        Findices[1].append(k)
                        Findices[2].append(m)
                    elif(tallyid == 2):
                        Cindices[0].append(j)
                        Cindices[1].append(k)
                        Cindices[2].append(m)			        
                elif(values[0].startswith("Flux")):
                    #bingo! We want the flux values
                    if(tallyid == 1):
                        Fval.append(float(values[1]))
                        Fvalsq.append(float(values[1]) * float(values[1]))
                    elif(tallyid == 2):
                        Cval.append(float(values[1]))
                        Cvalsq.append(float(values[1]) * float(values[1]))
                    i += 1 #Increment the line counter
            else: #NOT first file: just accumulate data!
                if(values[0].startswith("Flux")):
                    if(tallyid == 1):
                        Fval[i] += float(values[1])
                        Fvalsq[i] += float(values[1]) * float(values[1])
                    elif(tallyid == 2):
                        Cval[i] += float(values[1])
                        Cvalsq[i] += float(values[1]) * float(values[1])
                    i += 1 #Increment the line counter
            #print i
        #unset the first-file flag
        fileflag = 0
    else: #Uh-oh- this file seems to be missing!
        print "File " + `index` + "_tallies.out is missing!"
        realfiles -= 1 #decrement the number of 'real' files

Fcells = len(Findices[0])
Ccells = len(Cindices[0])

FRSD = [0]*Fcells
CRSD = [0]*Ccells

#Be sure to use REAL number of files to calc rsd:
#print "# of real files read: " + str(realfiles)
for i in range(Fcells):
	#Calc the RSD value 
	FRSD[i] = math.fabs( math.sqrt((realfiles/(realfiles - 1.0)) * math.fabs((Fvalsq[i])/realfiles - ((Fval[i])*(Fval[i]))/(realfiles*realfiles) ) ) / (Fval[i]/realfiles) )
	#print "RSD [" + str(g) + "][" + str(i) + "]: " + str(RSD[g][i])
	
for i in range(Ccells):
	#Calc the RSD value 
	CRSD[i] = math.fabs( math.sqrt((realfiles/(realfiles - 1.0)) * math.fabs((Cvalsq[i])/realfiles - ((Cval[i])*(Cval[i]))/(realfiles*realfiles) ) ) / (Cval[i]/realfiles) )
	#print "RSD [" + str(g) + "][" + str(i) + "]: " + str(RSD[g][i])

#All the error lines should be the same length
for i in range(Fcells):
	this_string = str(Findices[0][i]) + '\t' + str(Findices[1][i]) + '\t' + str(Findices[2][i]) + '\t' + str(FRSD[i]) #write cell indices, RSD
	dest_file1.write('\n' + this_string)

#All the error lines should be the same length
for i in range(Ccells):
	this_string = str(Cindices[0][i]) + '\t' + str(Cindices[1][i]) + '\t' + str(Cindices[2][i]) + '\t' + str(CRSD[i]) #write cell indices, RSD
	dest_file2.write('\n' + this_string)

dest_file1.close()
dest_file2.close() 

print "Done with real RSD calc.  See " + file_string1 + " and " + file_string2 + " for real RSD data."  




        
