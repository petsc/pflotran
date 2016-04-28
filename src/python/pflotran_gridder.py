#=========================================================
#
# Python script to create a structured mesh related files 
# to be read as an unstructured mesh in PFLOTRAN
# Author: Satish Karra, LANL
# Date: 07/08/2013
#
#=========================================================


import sys
import os
import numpy as np
import math

if len(sys.argv) < 10:
  sys.exit('ERROR: Usage - python gridder.py [No. of x verts][No. of y verts][No. of z verts][x_min][y_min][z_min][x_max][y_max][z_max]')
x_verts = int(sys.argv[1])
y_verts = int(sys.argv[2])
z_verts = int(sys.argv[3])
x_min = float(sys.argv[4])
y_min = float(sys.argv[5])
z_min = float(sys.argv[6])
x_max = float(sys.argv[7])
y_max = float(sys.argv[8])
z_max = float(sys.argv[9])


Total_verts = x_verts*y_verts*z_verts
N_cells = (x_verts-1)*(y_verts-1)*(z_verts-1)
delta_x = (x_max - x_min)/(x_verts - 1)
delta_y = (y_max - y_min)/(y_verts - 1)
delta_z = (z_max - z_min)/(z_verts - 1)

Coord = np.zeros((Total_verts,3),'float')
for i in range(1,x_verts+1):
    for j in range(1,y_verts+1):
        for k in range(1,z_verts+1):
            id = i+(j-1)*x_verts+(k-1)*x_verts*y_verts-1
            Coord[id,0] = (i-1)*delta_x
            Coord[id,1] = (j-1)*delta_y
            Coord[id,2] = (k-1)*delta_z
            
# Storing vertices in each element
# Assuming all elements are hexes
Vertices = np.zeros((N_cells,8),'int')

count = 0
for k in range(1,z_verts):
    for j in range(1,y_verts):
        for i in range(1,x_verts):
            id = i+(j-1)*x_verts+(k-1)*x_verts*y_verts
            Vertices[count,0] = id
            Vertices[count,1] = id+1
            Vertices[count,2] = id+x_verts+1
            Vertices[count,3] = id+x_verts
            Vertices[count,4] = id+x_verts*y_verts
            Vertices[count,5] = id+x_verts*y_verts+1
            Vertices[count,6] = id+x_verts*y_verts+x_verts+1
            Vertices[count,7] = id+x_verts*y_verts+x_verts
            count = count + 1
            
# Writing list of all vertices
print('--> Writing vertices')
fid = open('all.vset','w')
for i in range(1,Total_verts+1):
    fid.write('%i\n' %i)
fid.close()
print('--> Finished writing all.vset') 

# Writing mesh file
print('--> Writing usg file')
fid = open('usg.mesh','w')
fid.write('%i %i\n' %(N_cells,Total_verts))
for id in range(count):
    fid.write('H %i %i %i %i %i %i %i %i\n' %(Vertices[id,0],Vertices[id,1],Vertices[id,2],Vertices[id,3],Vertices[id,4],Vertices[id,5],Vertices[id,6],Vertices[id,7]))
for id in range(Total_verts):
    fid.write('%f %f %f\n' %(Coord[id,0],Coord[id,1],Coord[id,2]))
fid.close()
print('--> Finished writing usg.mesh')    

# Writing vertex numbers on faces
# Bottom (z=z_min)
print('--> Writing bottom vertices')
fid = open('bottom.vset','w')
for i in range(1,x_verts*y_verts+1):
    fid.write('%i\n' %i)
fid.close()
print('--> Finished writing bottom.vset') 


# Top (z=z_max)
print('--> Writing top vertices')
fid = open('top.vset','w')
for i in range(x_verts*y_verts*(z_verts-1),x_verts*y_verts*z_verts+1):
    fid.write('%i\n' %i)
fid.close()
print('--> Finished writing top.vset') 

# North (y=y_max)
print('--> Writing north vertices')
fid = open('north.vset','w')
for i in range(1,x_verts+1):
    for k in range(1,z_verts+1):
        j = y_verts
        id = i+(j-1)*x_verts+(k-1)*x_verts*y_verts
        fid.write('%i\n' %id)
fid.close()
print('--> Finished writing north.vset') 

# South (y=y_min)
print('--> Writing south vertices')
fid = open('south.vset','w')
for i in range(1,x_verts+1):
    for k in range(1,z_verts+1):
        j = 1
        id = i+(j-1)*x_verts+(k-1)*x_verts*y_verts
        fid.write('%i\n' %id)
fid.close()
print('--> Finished writing south.vset') 

# East (x=x_max)
print('--> Writing east vertices')
fid = open('east.vset','w')
for j in range(1,y_verts+1):
    for k in range(1,z_verts+1):
        i = x_verts
        id = i+(j-1)*x_verts+(k-1)*x_verts*y_verts
        fid.write('%i\n' %id)
fid.close()
print('--> Finished writing east.vset') 

# West (x=x_min)
print('--> Writing west vertices')
fid = open('west.vset','w')
for j in range(1,y_verts+1):
    for k in range(1,z_verts+1):
        i = 1
        id = i+(j-1)*x_verts+(k-1)*x_verts*y_verts
        fid.write('%i\n' %id)
fid.close()
print('--> Finished writing west.vset') 

# Create a subdirectory
d = 'dat'
import shutil
if os.path.isdir(d): #check if d exists
    shutil.rmtree(d) #remove old directory
os.mkdir(d)          #create new directory

# Move *.vset *.mesh files to new subdirectory
cmd = 'mv' + ' *.vset' + ' usg.mesh' + ' %s/.'%d
failure = os.system(cmd)
if failure:
    print 'Unable to move *.vset, *.mesh files to subdirectory'; sys.exit(1)
print('--> Finished with moving files to dat directory')



 
