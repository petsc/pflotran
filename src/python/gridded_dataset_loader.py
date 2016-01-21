import sys
from h5py import *
import numpy

filename = 'dataset.h5'
h5file = File(filename,mode='w')

# 2d surface
h5grp = h5file.create_group('test_surface')

nx = 3
ny = 3
nz = 1
nt = 50

h5grp.attrs['Dimension'] = numpy.string_('XY')
h5grp.attrs['Discretization'] = [2.5,2.5]
h5grp.attrs['Origin'] = [0.,0.]
h5grp.attrs['Max Buffer Size'] = [4]

rarray = numpy.zeros((nt),'=f8')
for t in range(nt):
  rarray[t] = 10.*float(t)
h5dset = h5grp.create_dataset('Times', data=rarray)

rarray = numpy.zeros((nx,ny,nt),'=f8')

for t in range(nt):
  for j in range(ny):
    for i in range(nx):
      rarray[i][j][t] = float(i) + 1.e-2*float(j) + 2. + \
                        1.e2*float(t)
h5dset = h5grp.create_dataset('Data', data=rarray)

# 1d line
h5grp = h5file.create_group('x_line')

nx = 5

h5grp.attrs['Dimension'] = numpy.string_('X')
h5grp.attrs['Discretization'] = [2.]
h5grp.attrs['Origin'] = [-1.]
h5grp.attrs['Max Buffer Size'] = [2]

rarray = numpy.zeros(nx,'=f8')

for i in range(nx):
  rarray[i] = 5. - 0.1*float(i)
h5dset = h5grp.create_dataset('Data', data=rarray)
 
# 1d line
h5grp = h5file.create_group('y_line')

nx = 1
ny = 5
nz = 1
nt = 30

h5grp.attrs['Dimension'] = numpy.string_('Y')
h5grp.attrs['Discretization'] = [2.]
h5grp.attrs['Origin'] = [-1.] 
h5grp.attrs['Transient'] = [True]
h5grp.attrs['Max Buffer Size'] = [2]

rarray = numpy.zeros((nt),'=f8')
for t in range(nt):
  rarray[t] = 5.*float(t)
h5dset = h5grp.create_dataset('Times', data=rarray)

rarray = numpy.zeros((ny,nt),'=f8')

for t in range(nt):
  for j in range(ny):
    rarray[j][t] = 5. + 0.1*float(j) + float(t)
h5dset = h5grp.create_dataset('Data', data=rarray)
 
h5file.close()

print('done')
