import sys
import os
import math
import numpy as npy

try:
    #python 3-friendly version
    exec(open("parameters.inp").read(), globals())
except:
    print( "Uh-oh. Where is parameters.inp?")
    exit(1)

#Now that we have the content of parameters, get the number of files
#(same as no. of different coarse grids)
n_files = len(cg_set)
n_relfs = len(rf_set)

#Figure out how many total cases (xsect combos) we should have for each grid/relax factor combo:
n_cases = len(st1_set) * len(st2_set) * len(c1_set) * len(ss12_c1_set) * len(c2_set)
print( str(n_cases) + " cases found in parameters.inp!\n" )

#There will be n_relfs rows, and about a zillion cases in each row
#We'll just keep appending all of the grid-size file data onto the end of the rows
specrads = [[] for b in range(n_relfs)]

for a in range(n_files):
    file_name = "Implicit_2g_delta_" + str(cg_set[a]) + ".out"
    try:
        this_file = open(file_name, 'r')
    except:
        print "File " + file_name + " not found! :("
        exit(1)
    print( "Processing file " + file_name )
    index = 0
    for line in this_file:
        #strip the newline and any extra righthand spaces
        line = line.rstrip(' \n')
        #split by spaces:
        entries = line.split(' ')
        if(entries[0].startswith("#")):
            #print( "We don't want this one: " + entries[0] )
            continue
        else:
            #Does the number of entries match the number of cases (+1 for the rel. factor)?
            if((len(entries)-1) != n_cases):
                print("parameters.inp suggests " + str(n_cases) + "cases; only " + str(len(entries)-1) + "found!")
                exit(1)
            #print( "We DO want this one: " + str(len(entries)))
            #Stash everything but the FIRST entry (value of relf)
            for j in range(len(entries)-1):
                specrads[index].append(entries[j+1])
        index=index+1
    this_file.close()

#Once we have ALL the data, we need to search through each COLUMN and save two
#values: (1) the "optimum" relaxation factor (the value at which the spectral 
#radius is lowest), and (2) the spectral radius value for relf=1.0.  
final_coords = [[] for i in range(2)]
#final_coords = npy.zeros((2, n_cases*n_files))

#for all of the columns...
#for k in range(len(specrads[0])):
for k in range(n_cases*n_files):
    #for all of the rows...
    opt_relf = 0.0
    min_spec = 1e25
    for l in range(n_relfs):
        #Is this spectral radius the new min?
        if(float(specrads[l][k]) < float(min_spec)):
            #If yes, store it!
            #print("values: " + str(min_spec) + " " + str(specrads[l][k]) + " " + str(rf_set[l]))
            min_spec = specrads[l][k]
            opt_relf = rf_set[l]
    #Cases where RF=1.0 specrad is > 1.2 are basically garbage :D
    if(float(specrads[n_relfs-1][k]) < 1.2):
        #Append this coordinate pair! (last row of specrads will always be RF=1.0):         
        final_coords[0].append(specrads[n_relfs-1][k])
        final_coords[1].append(opt_relf)

print len(final_coords[1])
print len(final_coords[0])
#Use final coordinate data to build empirical fits (linear, quad):
final_coord0 = npy.asarray(final_coords[0], dtype=npy.float32)
final_coord1 = npy.asarray(final_coords[1], dtype=npy.float32)
lin_coeffs, lin_resid, lin_rank, lin_sing, lin_cond = npy.polyfit(final_coord0, final_coord1, 1, full=True)

print lin_coeffs
print lin_resid
#Now, dump it all into a file!
outfile = open("Opt_relax_data_fit.out", 'w')

print "Printing results of " + str(len(final_coords[0])) + " FA simulations..."

outfile.write("#Spec Rad, w=1 | Optimal w\n")
for i in range(len(final_coords[0])):
    outfile.write(str(final_coords[0][i]) + '\t' + str(final_coords[1][i]) + '\n')          
outfile.close()

print "                                               ...Finished!!"      
