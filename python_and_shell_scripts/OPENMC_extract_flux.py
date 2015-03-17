from __future__ import division
import sys
import os.path
import math

if(len(sys.argv) < 3):
    print( "Usage: python OPENMC_extract_flux.py <TALLY> <OUTPUT NAME>")
    sys.exit()

tallyid = int(sys.argv[1]) #Tally id to look for-- passed from terminal
file_string1 = sys.argv[2] #Name of output to be written to

print( "Writing data from tally " + sys.argv[1] + " to file " + sys.argv[2] + ". ")


index = 0
printflag = 1
fileflag = 1
Findices = [[],[],[]]
Gindices = []
ROIindices = [[],[],[]]

Fval = []
Ferror = []

flag = 0

i = 0
if os.path.isfile("tallies.out") is True:
    print( "Found file tallies.out!")
    this_file = open("tallies.out", 'r')
    for line in this_file:
        #split into values by whitespace
        values = line.split()
        if(len(values) == 0): 
            continue
        elif(values[0].startswith("===")):
            #This is a tally identifier line!
            if(int(values[2].strip(':')) == tallyid):
                flag = 1
                continue   
            else:
                flag = 0   #This is not a tally we want- set reading flag to 0
                continue 
        elif(values[0].startswith("Total") or values[0].startswith("Vol")):
            continue    #skip this line
        elif(values[0].startswith("Incoming") and values[1].startswith("Energy")): #There's an energy bin on this line- take average of bin bounds to get center
            g_en = 0.5 * (float(values[2].strip('[,)')) + float(values[3].strip('[,)')))
            if(flag == 1):
                Gindices.append(g_en)
        elif(values[0].startswith("Mesh")):
            #strip all commas, parens. Grab indices!
            if(flag == 1):
                j = int(values[2].strip('(,)'))
                k = int(values[3].strip('(,)'))
                if(len(values) > 4):
                    m = int(values[4].strip('(,)'))
                else:
                    m = 1
                #store indices:
                Findices[0].append(j)
                Findices[1].append(k)
                Findices[2].append(m)
            else:
                continue				        
        elif(values[0].startswith("Flux")):
            #bingo! We want the flux values
            if(flag == 1):
                if(printflag==1):
                    print( "It's a flux tally!")
                    printflag = 0
                Fval.append(float(values[1]))
                Ferror.append(float(values[3]))
            i += 1 #Increment the line counter
        elif(values[0].startswith("Fission")):
            if(flag == 1):
                if(printflag==1):
                    print( "It's a fission tally!")
                    printflag = 0
                Fval.append(float(values[2]))
                Ferror.append(float(values[4]))
            i += 1 #Increment the line counter


Fcells = len(Findices[0])
n_grp = int(len(Gindices)/float(Fcells))
if(n_grp == 0):
    #Energy-integrated, just set to 1
    n_grp = 1
print( "Total number of groups: " + str(n_grp))

if(Fcells > 0):
    dest_file1 = open(file_string1, 'w')
    title_string = "#X\tY\tZ\t"
    for a in range(n_grp):
        title_string += "G" + str(a+1) + " flux \t+/- Std. Dev.\t\t" 
    dest_file1.write(title_string)
    print( "Created output file " + file_string1)
    
    #All the error lines should be the same length
    for i in range(Fcells):
	    this_string = str(Findices[0][i]) + '\t' + str(Findices[1][i]) + '\t' + str(Findices[2][i]) + '\t' 
	    for j in range(n_grp):
	        this_string = this_string  + ("%7.6e") % (Fval[i*n_grp+j]) + "\t" + ("%7.6e") % (Ferror[i*n_grp+j]) + '\t' #append cell indices, flux, app rsd for this group
	    dest_file1.write('\n' + this_string)

    dest_file1.close()

if(Fcells > 0):
    print( "Done printing tally " + sys.argv[1] + "!  See " + file_string1 + " for output." )
else:
    print( "Tally " + sys.argv[1] + " not found in tallies.out!")
     





        
