import sys
from h5py import *
import numpy
import random

def write_array(fid,array):
  num_per_line = 4
  array_size = len(array)
  for i in range(array_size):
    if i % num_per_line == 0:
      fid.write('  ')
    fid.write('%.8e ' % array[i])
    if (i+1) % num_per_line == 0 and i+1 != array_size:
        fid.write('\\\n')
  fid.write('\n')

num_times = 21
num_values_per_time = 11
tdata = [0.]*num_times
pdata = [0.]*num_times*num_values_per_time
sdata = [0.]*num_times*num_values_per_time

max_pressure = [6.e5]*num_times
min_pressure = [101325.]*num_times
max_porosity = [0.3]*num_times
min_porosity = [0.05]*num_times
max_pressure[0] = 2.e5
max_pressure[1] = 2.e5
max_porosity[1] = 0.35
for i in range(17,num_times):
  min_porosity[i] = 0.15

for itime in range(num_times):
  tdata[itime] = float(itime)/(num_times-1)
  rand = [0.]*num_values_per_time
  sum = 0.
  for j in range(1,num_values_per_time):
    rand[j] = random.random()
    sum += rand[j]
  pres_diff = max_pressure[itime]-min_pressure[itime]
  pres_scale = pres_diff/sum
  por_diff = max_porosity[itime]-min_porosity[itime]
  for j in range(num_values_per_time):
    index = itime*num_values_per_time+j
    if j == 0:
      pdata[index] = rand[j]*pres_scale + min_pressure[itime]
    else:
      pdata[index] = rand[j]*pres_scale + pdata[index-1]
    sdata[index] = (pdata[index]-min_pressure[itime])/pres_diff* \
                   por_diff + min_porosity[itime]

filename = 'pflotran_closure_11x11.dat'
file = open(filename,'w')
file.write('NUM_TIMES %d\n' % num_times)
file.write('NUM_VALUES_PER_TIME %d\n' % num_values_per_time)
file.write('TIME_UNITS y\n')
file.write('TIME \\\n')
write_array(file,tdata)
file.write('PRESSURE \\\n')
write_array(file,pdata)
file.write('POROSITY \\\n')
write_array(file,sdata)
file.close()

print('done')
