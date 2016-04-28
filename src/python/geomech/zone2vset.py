#--------------------------------------------------
# Program to read boundary node ids from zone file
# Then write to vset file
# Author: Satish Karra
# Date: July 11, 2013
# Email: satkarra@lanl.gov
#--------------------------------------------------

import numpy as np
import sys

if len(sys.argv) < 3:
    sys.exit('ERROR: Usage - python zone2vset.py [zone_filename][vset_filename]')
zone_file = sys.argv[1]
vset_file = sys.argv[2]

# Opening zone file
print ('--> Reading zone file')

f = open(zone_file,'r')
f.readline()
f.readline()
f.readline()

NumNodes = int(f.readline())
Node_array = np.zeros(NumNodes,'int')

if (NumNodes < 10):
	g = f.readline()
	node_array = g.split()
# Convert string to integer array
	node_array = [int(id) for id in node_array]
	Node_array = node_array
else:
	for i in range(NumNodes/10 + 1):
		g = f.readline()
		node_array = g.split()
# Convert string to integer array
		node_array = [int(id) for id in node_array]
		if (NumNodes-10*i < 10):
			for j in range(NumNodes%10):
				Node_array[i*10 + j] = node_array[j]
		else:
			for j in range(10):
				Node_array[i*10 + j] = node_array[j]

print('--> Finished with zone file')



# Write to vset file
print('--> Writing to vset file')
fid = open(vset_file,'w')
for i in range(NumNodes):
    fid.write('%i \n'  %Node_array[i])

print('--> Finished writing to vset file')
