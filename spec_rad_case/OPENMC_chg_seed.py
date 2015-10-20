import sys
import os
import math

seed = int(sys.argv[1]) #Seed value to insert-- passed from bash script

edit_file = open('settings.xml', 'r')
temp_file = open('temp', 'w')

for line in edit_file:
    stripped = line.strip()   #kill spaces
    if(stripped.startswith("<seed>")):
        temp_file.write("    <seed>" + str(seed) + "</seed>\n")
    else:
        temp_file.write(line)
        
        
edit_file.close()
temp_file.close()

os.rename('temp', 'settings.xml')
