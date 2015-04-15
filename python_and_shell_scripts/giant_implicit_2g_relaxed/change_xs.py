import sys
import os
import math

cg_index = int(sys.argv[1])-1 #coarse grid index-- passed from bash script
rf_index = int(sys.argv[2])-1 #relax fact index-- passed from bash script
st1_index = int(sys.argv[3])-1 #sigmaT1 index-- passed from bash script
st2_index = int(sys.argv[4])-1 #sigmaT2 index-- passed from bash script
c1_index = int(sys.argv[5])-1 #sigmaS11 index-- passed from bash script
ss12_c1_index = int(sys.argv[6])-1 #sigmaS12 index-- passed from bash script
c2_index = int(sys.argv[7])-1 #sigmaS22 index-- passed from bash script

edit_file = open('2g_test_input.inp', 'r')
temp_file = open('temp', 'w')
try:
    execfile("parameters.inp")
except:
    print( "Uh-oh. Where is parameters.inp?")

this_cg = cg_set[cg_index]
this_rf = rf_set[rf_index]
this_st1 = st1_set[st1_index]
this_st2 = st2_set[st2_index]
this_c1 = c1_set[c1_index]
this_ss12_c1 = ss12_c1_set[ss12_c1_index]
this_c2 = c2_set[c2_index]

for line in edit_file:
    stripped = line.strip()   #kill spaces
    if(stripped.startswith("sigT = ")): #Write total xsects!
        temp_file.write("sigT = [" + str(this_st1) + ", " + str(this_st2) + "]\n")
    elif(stripped.startswith("h_j = ")): #Write coarse grid!
        temp_file.write("h_j = " + str(this_cg) + "\n")
    elif(stripped.startswith("sigS = ")): #Write scattering xsects!
        ss12 = this_ss12_c1 * this_c1 * this_st1
        ss11 = this_c1 * this_st1 - ss12
        ss22 = this_st2 * this_c2
        temp_file.write("sigS = [[" + str(ss11) + ", " + str(ss12) + "],[0.0, " + str(ss22) + "]]\n")
    elif(stripped.startswith("rel_fact = ")): #Write relaxation factor!
        temp_file.write("rel_fact = " + str(this_rf) + "\n")
    else:
        temp_file.write(line)
        
        
edit_file.close()
temp_file.close()

os.rename('temp', '2g_test_input.inp')
