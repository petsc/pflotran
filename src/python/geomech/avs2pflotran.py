#-------------------------------------------------
# Program to read avs file and convert to PFLOTRAN
# implicit grid format
# Author: Satish Karra
# Email: satkarra@lanl.gov
# Date: July 11, 2013
#-------------------------------------------------

import numpy as np
import sys

if len(sys.argv) < 2:
    sys.exit('Error: Usage - python avs2pflotran [avs_filename]')
avs_file = sys.argv[1]

# Opening avs file
print('--> Opening avs file')
f = open(avs_file,'r')
g = f.readline()
g = g.split()

Num_Nodes = int(g.pop(0))
Num_Elem = int(g.pop(0))


Vertices = np.zeros((Num_Nodes,3),'float')
for i in range(Num_Nodes):
    g = f.readline()
    g = g.split()
    Vertices[i,0] = float(g.pop(1))
    Vertices[i,1] = float(g.pop(1))
    Vertices[i,2] = float(g.pop(1))

# elem_type
# tet = 4
# brick = 8


Elements = np.zeros((Num_Elem,9),'int')
for i in range(Num_Elem):
    g = f.readline()
    g = g.split()
    elem_type = g.pop(2)
    if (elem_type == 'tet'):
        Elements[i,0] = 4
    for j in range(Elements[i,0]):
        Elements[i,j+1] = g.pop(2)

print('--> Done with avs file')

# Writing to pflotran implicit grid format
# Note that vertex 4 is written first and vertex 1 is written last
# This is to ensure the right hand rule for tetrahedron is preserved
# Check with Carl why lagrit prints tets in clockwise direction
print('--> Writing pflotran usg file')
fid = open('usg.mesh','w')
fid.write('%i %i\n' %(Num_Elem,Num_Nodes))
for i in range(Num_Elem):
    if (Elements[i,0] == 4):
        fid.write('T %i %i %i %i\n' %(Elements[i,4],Elements[i,2],Elements[i,3],Elements[i,1]))
    else:
        sys.exit('Error: Invalid element type')

# Writing coordinates
for i in range(Num_Nodes):
    fid.write('%f %f %f\n' %(Vertices[i,0],Vertices[i,1],Vertices[i,2]))
fid.close()
print('--> Finished writing usg.mesh')


# Writing all vertices
print('--> Writing vertices')
fid = open('all.vset','w')
for i in range(1,Num_Nodes+1):
    fid.write('%i \n' %i)
fid.close()
print('--> Finished writing all.vset')


