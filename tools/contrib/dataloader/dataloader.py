# Author: Glenn Hammond
# Date: 09/06/13
# Copyright: Pacific Northwest National Laboratory

import sys
from h5py import *
import numpy

meters_per_day_to_meters_sq = 0.89e-3/997.32/9.81/24./3600.

filename = 'permeability.h5'
h5file = File(filename,mode='w')

nx = 3
ny = 3
nz = 3

print 'allocating cell index array....'
iarray = numpy.zeros((nx*ny*nz),'=i4')

print 'setting cell indices....'
# add cell ids to file
for i in range(nx*ny*nz):
  iarray[i] = i+1
dataset_name = 'Cell Ids'
h5dset = h5file.create_dataset(dataset_name, data=iarray)

print 'allocating data array....'
a = numpy.zeros((nx*ny*nz),'=f8')

filename = 'dataloader.out'
fout = open(filename,'w')
print 'reading data....'
filename = 'perm_test.dat'
f = open(filename,'r')
count = 0
min = 1e20
max = -1.e20
ave = 0.
# The ordering of these loops is for first k, then j, then i, i.e.
#   index = k + j*nz + i*ny*nz
# which is how the test file is ordered
for i in range(nx):
  for j in range(ny):
    for k in range(nz):
      # this index will convert to i, j, k format
      index = i + j*nx + k*nx*ny
      s = f.readline()
#### COMMENT OUT LINES BETWEEN IF I,J,K INDICES NOT INCLUDED IN TEXT INPUT FILE
      # since the string has coordinates at the beginning, split it based on
      # white space
      w = s.split() # white space is the default 
      # set string equal to second substring
      s = w[1]
#### COMMENT OUT LINES BETWEEN IF I,J,K INDICES NOT INCLUDED IN TEXT INPUT FILE
      value = float(s)
      if value < min:
        min = value
      if value > max:
        max = value
      ave += value
      a[index] = value* meters_per_day_to_meters_sq 
      count += 1
      if count % 10000000 == 0:
        print count
dataset_name = 'Permeability'
h5dset = h5file.create_dataset(dataset_name, data=a)
print 'min: %e' % min
print 'max: %e' % max
print 'ave: %e' % (ave/float(count))
print '\n'
fout.write('min: %e\n' % min)
fout.write('max: %e\n' % max)
fout.write('ave: %e\n' % (ave/float(count)))
fout.write('\n')
print 'done with ', dataset_name
f.close()
fout.close()
  
h5file.close()
print 'done with everything'
  
