from __future__ import division
import sys
import os.path
import math

file_string1 = "MC_flux.out"

dest_file1 = open(file_string1, 'w')

print "Created output file " + file_string1
dest_file1.write("Mesh indices, flux, apparent RSD: ")

index = 0
fileflag = 1
Findices = [[],[],[]]

Fval = []
Ferror = []

tallyid = 0

i = 0
if os.path.isfile("tallies.out") is True:
    print "Found file tallies.out!"
    this_file = open("tallies.out", 'r')
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
            else:
                break   #This is not a tally we want!
        elif(values[0].startswith("Total") or values[0].startswith("Vol")):
            continue    #skip this line
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
        elif(values[0].startswith("Flux")):
            #bingo! We want the flux values
            if(tallyid == 1):
                Fval.append(float(values[1]))
                Ferror.append(float(values[3]))
            i += 1 #Increment the line counter


Fcells = len(Findices[0])

FRSD = [0]*Fcells

#All the error lines should be the same length
for i in range(Fcells):
	this_string = str(Findices[0][i]) + '\t' + str(Findices[1][i]) + '\t' + str(Findices[2][i]) + '\t' + str(Fval[i]) + '\t' + str(Ferror[i]/Fval[i]) #write cell indices, flux, app rsd
	dest_file1.write('\n' + this_string)

dest_file1.close()

print "Done printing MC flux!  See " + file_string1 + " for output."  




        
