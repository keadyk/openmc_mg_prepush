import sys
import os
import math

this_cg_index = int(sys.argv[1])-1 #coarse grid index-- passed from bash script

try:
    execfile("parameters.inp")
except:
    print "Uh-oh. Where is parameters.inp?"

try:
    execfile("2g_test_input.inp")
except:
    print "Uh-oh. Where is 2g_test_input.inp?"

file_name = "Implicit_2g_delta_" + str(cg_set[this_cg_index]) + ".out"

grid_file = open(file_name, 'w')

#Write the header line with the stuff that doesn't change:
grid_file.write("#Slab size " + str(X) + " cm, M=" + str(M) + ", h_k=" + str(h_k) + "\n")
grid_file.write("#")
#Write all the case infos
for a in range(len(st1_set)):
    for b in range(len(st2_set)):
        for c in range(len(c1_set)):
            for d in range(len(ss12_ss1_set)):
                for e in range(len(c2_set)): 
                    ss12 = ss12_ss1_set[d] * c1_set[c] * st1_set[a]
                    ss11 = c1_set[c]*st1_set[a] - ss12
                    ss22 = st2_set[b] * c2_set[e]
                    grid_file.write(" " + str(st1_set[a]) + " " + str(st2_set[b]) + " " + str(ss11) + " " + str(ss12) + " " + str(ss22) + " |")
grid_file.write("\n")
grid_file.close()
