import sys
import math
from h5py import *
import numpy

filename = 'parameters.h5'
#filename = 'input_small_100.h5'
h5file = File(filename,mode='w')

nx = 170
ny = 200
nz = 40
n = nx*ny*nz

#meters_per_day_to_meters_sq = 0.89e-3/997.32/9.81/24./3600.
#meters_per_day_to_meters_sq = 1.

iarray = numpy.zeros((nx*ny*nz),'=i4')

# add cell ids to file
for i in range(n):
  iarray[i] = i+1
dataset_name = 'Cell Ids'
h5dset = h5file.create_dataset(dataset_name, data=iarray)

a = numpy.zeros(n,'=f8')

filename = 'permeability.final'
f = open(filename,'r')
# remove top 3 lines
s = f.readline()
s = f.readline()
s = f.readline()

ln_average = 0.
average = 0.
min = 1e20
max = -1e20
for i in range(n):
  s = f.readline()
  w = s.split()
  ln_value = float(w[0])
  value = math.exp(ln_value)
  if value < min:
    min = value
  if value > max: 
    max = value
  a[i] = value
  average += value
  ln_average += ln_value
dataset_name = 'Permeability'
h5dset = h5file.create_dataset(dataset_name, data=a)

average /= n
ln_average /= n

stdev = 0.
ln_stdev = 0.
for i in range(n):
  stdev += (a[i]-average)*(a[i]-average)
  value = math.log(a[i])
  ln_stdev += (value-ln_average)*(value-ln_average)
stdev /= n
ln_stdev /= n

stdev = math.sqrt(stdev)
ln_stdev = math.sqrt(ln_stdev)

print '\nPermeability Statistics'
print 'average: %e  stdev: %e' % (average,stdev)
print 'ln_average: %e  ln_stdev: %e' % (ln_average,ln_stdev)
print 'min: %e' % min
print 'max: %e' % max

print '\ndone with ', dataset_name
f.close()
  
filename = 'porosity.final'
f = open(filename,'r')
# remove top 3 lines
s = f.readline()
s = f.readline()
s = f.readline()

num_below_05 = 0
num_below_0 = 0
average = 0.
min = 1e20
max = -1e20
for i in range(n):
  s = f.readline()
  w = s.split()
  value = float(w[0])
  if value < 0.05:
    if value < 0.:
      num_below_0 += 1
    num_below_05 += 1
 #   print 'porosity value %f truncated to 0.05' % value
    value = 0.05
  if value < min:
    min = value
  if value > max: 
    max = value
  a[i] = value
  average += value
dataset_name = 'Porosity'
h5dset = h5file.create_dataset(dataset_name, data=a)

average /= n

stdev = 0.
for i in range(n):
  stdev += (a[i]-average)*(a[i]-average)
stdev /= n
stdev = math.sqrt(stdev)

print '\nPorosity Statistics'
print 'average: %e  :stdev %e' % (average,stdev)
print 'min: %e' % min
print 'max: %e' % max
print 'Number of values below 0.05: %d' % num_below_05
print 'Number of negative values: %d' % num_below_0

print '\ndone with ', dataset_name
f.close()
  
  
h5file.close()
print 'done with everything'
  
