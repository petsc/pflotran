#!/usr/bin/python

#------------------------------------------------------------------------
# Program to generate unstructured explicit format for a structured mesh
# for PFLOTRAN (www.pflotran.org)
# Author: Satish Karra, LANL
# Date: Sept. 2, 2013
#------------------------------------------------------------------------

import sys
import os
import numpy as np
import math

# Reading the data
if len(sys.argv) < 10:
  sys.exit('ERROR: Usage - python structured_to_explicit.py [No. of x verts]\
  	        [No. of y verts][No. of z verts][x_min][y_min][z_min][x_max][y_max][z_max]')
x_verts = int(sys.argv[1])
y_verts = int(sys.argv[2])
z_verts = int(sys.argv[3])
x_min = float(sys.argv[4])
y_min = float(sys.argv[5])
z_min = float(sys.argv[6])
x_max = float(sys.argv[7])
y_max = float(sys.argv[8])
z_max = float(sys.argv[9])

# Num of cells, increments in x, y, z
N_cells = (x_verts-1)*(y_verts-1)*(z_verts-1)
delta_x = (x_max - x_min)/(x_verts - 1)
delta_y = (y_max - y_min)/(y_verts - 1)
delta_z = (z_max - z_min)/(z_verts - 1)

vol = delta_x*delta_y*delta_z

count = 0
Coord = np.zeros((N_cells,3),'float')
for k in range(1,z_verts):
    for j in range(1,y_verts):
        for i in range(1,x_verts):
            Coord[count,0] = x_min + delta_x/2.0 + (i-1)*delta_x
            Coord[count,1] = y_min + delta_y/2.0 + (j-1)*delta_y
            Coord[count,2] = z_min + delta_z/2.0 + (k-1)*delta_z
            count = count + 1


x_conn = (x_verts - 2)*(y_verts - 1)*(z_verts - 1)
y_conn = (y_verts - 2)*(x_verts - 1)*(z_verts - 1)
z_conn = (z_verts - 2)*(y_verts - 1)*(x_verts - 1)

Num_conn = x_conn + y_conn + z_conn

x_area = delta_z*delta_y
y_area = delta_x*delta_z
z_area = delta_x*delta_y

count = 0
Conn_coord = np.zeros((Num_conn,3),'float')
ids = np.zeros((Num_conn,2),'int')

for k in range(1,z_verts):
    for j in range(1,y_verts):
        for i in range(1,x_verts-1):
            ids[count,0] = i + (j-1)*(x_verts-1) + (k-1)*(x_verts-1)*(y_verts-1)
            ids[count,1] = ids[count,0] + 1
            count = count + 1

for k in range(1,z_verts):
    for i in range(1,x_verts):
        for j in range(1,y_verts-1):
            ids[count,0] = i + (j-1)*(x_verts-1) + (k-1)*(x_verts-1)*(y_verts-1)
            ids[count,1] = ids[count,0] + (x_verts-1) 
            count = count + 1

for i in range(1,x_verts):
    for j in range(1,y_verts):
        for k in range(1,z_verts-1):
            ids[count,0] = i + (j-1)*(x_verts-1) + (k-1)*(x_verts-1)*(y_verts-1)
            ids[count,1] = ids[count,0] + (x_verts-1)*(y_verts-1) 
            count = count + 1

count = 0
for i in range(Num_conn):
    Conn_coord[i,0] = 0.5*(Coord[ids[i,0]-1,0] + Coord[ids[i,1]-1,0])
    Conn_coord[i,1] = 0.5*(Coord[ids[i,0]-1,1] + Coord[ids[i,1]-1,1])
    Conn_coord[i,2] = 0.5*(Coord[ids[i,0]-1,2] + Coord[ids[i,1]-1,2])

# Areas
areas = np.zeros((Num_conn,1),'float')
for i in range(Num_conn):
    a = np.array([Coord[ids[i,0]-1,0],Coord[ids[i,0]-1,1],Coord[ids[i,0]-1,2]])
    b = np.array([Coord[ids[i,1]-1,0],Coord[ids[i,1]-1,1],Coord[ids[i,1]-1,2]])
    normal = b-a
    normal = normal/np.sqrt(np.inner(normal,normal))
    if (normal[1] == 1.0):
        areas[i] = x_area
    elif (normal[1] == 1.0):
        areas[i] = y_area
    else:
        areas[i] = z_area

#------------------------------------------------------------------------
# Writing explicit file
print('--> Writing uge file')
fid = open('grid.uge','w')
fid.write('%s  %i\n' %('CELLS',N_cells))
for id in range(N_cells):
    fid.write('%i %f %f %f %f \n' %(id+1, Coord[id,0],Coord[id,1],Coord[id,2],vol))
fid.write('%s  %i\n' %('CONNECTIONS',Num_conn))
for id in range(Num_conn):
    fid.write('%i %i %f %f %f %f\n' %(ids[id,0],ids[id,1],\
    	Conn_coord[id,0],Conn_coord[id,1],Conn_coord[id,2],areas[id]))
fid.close()
print('--> Finished writing grid.uge')   

#------------------------------------------------------------------------
# Boundary connections
# West connections x = xmin
id = np.zeros(((y_verts-1)*(z_verts-1),1),'int')
for i in range(1,(y_verts-1)*(z_verts-1)+1):
    id[i-1] = (i-1)*(x_verts-1) + 1

print('--> Writing west.ex file')
fid = open('west.ex','w')
fid.write('%s  %i\n' %('CONNECTIONS',(y_verts-1)*(z_verts-1)))
for i in range((y_verts-1)*(z_verts-1)):
    fid.write('%i %f %f %f %f\n' %(id[i],x_min,Coord[id[i]-1,1],\
    	Coord[id[i]-1,2],x_area))
fid.close()
print('--> Finished writing west.ex')    

#------------------------------------------------------------------------
# West connections x = x_max
id = np.zeros(((y_verts-1)*(z_verts-1),1),'int')
for i in range(1,(y_verts-1)*(z_verts-1)+1):
    id[i-1] = i*(x_verts-1)

print('--> Writing east.ex file')
fid = open('east.ex','w')
fid.write('%s  %i\n' %('CONNECTIONS',(y_verts-1)*(z_verts-1)))
for i in range((y_verts-1)*(z_verts-1)):
    fid.write('%i %f %f %f %f\n' %(id[i],x_max,Coord[id[i]-1,1],\
    	Coord[id[i]-1,2],x_area))
fid.close()
print('--> Finished writing east.ex')    

#------------------------------------------------------------------------
# south connections y = ymin
id = np.zeros(((x_verts-1)*(z_verts-1),1),'int')
count = 0
for k in range(1,z_verts):
    for i in range(1,x_verts):
        id[count] = i + (k-1)*(y_verts-1)*(x_verts-1)
        count = count + 1

print('--> Writing south.ex file')
fid = open('south.ex','w')
fid.write('%s  %i\n' %('CONNECTIONS',(x_verts-1)*(z_verts-1)))
for i in range((x_verts-1)*(z_verts-1)):
    fid.write('%i %f %f %f %f\n' %(id[i],Coord[id[i]-1,0],y_min,\
    	Coord[id[i]-1,2],y_area))
fid.close()
print('--> Finished writing south.ex')  

#------------------------------------------------------------------------
# north connections y = ymax
id = np.zeros(((x_verts-1)*(z_verts-1),1),'int')
count = 0
for k in range(1,z_verts):
    for i in range(1,x_verts):
        id[count] = i + (k-1)*(y_verts-1)*(x_verts-1) + (x_verts-1)*(y_verts-2)
        count = count + 1

print('--> Writing north.ex file')
fid = open('north.ex','w')
fid.write('%s  %i\n' %('CONNECTIONS',(x_verts-1)*(z_verts-1)))
for i in range((x_verts-1)*(z_verts-1)):
    fid.write('%i %f %f %f %f\n' %(id[i],Coord[id[i]-1,0],y_max,\
    	Coord[id[i]-1,2],y_area))
fid.close()
print('--> Finished writing north.ex')  

#------------------------------------------------------------------------
# bottom connections z = zmin
id = np.zeros(((x_verts-1)*(y_verts-1),1),'int')
count = 0
for j in range(1,y_verts):
    for i in range(1,x_verts):
        id[count] = i + (j-1)*(x_verts-1)
        count = count + 1

print('--> Writing bottom.ex file')
fid = open('bottom.ex','w')
fid.write('%s  %i\n' %('CONNECTIONS',(x_verts-1)*(z_verts-1)))
for i in range((x_verts-1)*(y_verts-1)):
    fid.write('%i %f %f %f %f\n' %(id[i],Coord[id[i]-1,0],\
    	Coord[id[i]-1,1],z_min,z_area))
fid.close()
print('--> Finished writing bottom.ex')  

#------------------------------------------------------------------------
# top connections z = zmax
id = np.zeros(((x_verts-1)*(y_verts-1),1),'int')
count = 0
for j in range(1,y_verts):
    for i in range(1,x_verts):
        id[count] = i + (j-1)*(x_verts-1) + (x_verts-1)*(y_verts-1)*(z_verts-2)
        count = count + 1

print('--> Writing top.ex file')
fid = open('top.ex','w')
fid.write('%s  %i\n' %('CONNECTIONS',(x_verts-1)*(z_verts-1)))
for i in range((x_verts-1)*(y_verts-1)):
    fid.write('%i %f %f %f %f\n' %(id[i],Coord[id[i]-1,0],\
    	Coord[id[i]-1,1],z_max,z_area))
fid.close()
print('--> Finished writing top.ex')  


