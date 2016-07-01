import sys
from h5py import *
import numpy as np

filename = 'dataset.h5'
h5file = File(filename,mode='w')

L = 1.

# 1d line in x
# Temperature boundary condition; SOUTH
h5grp = h5file.create_group('x_line_node_centered_south')
h5grp.attrs['Dimension'] = np.string_('X')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [2.]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
nx = 2
rarray = np.zeros(nx,'=f8')
for i in range(nx):
  # T = x
  x = (float(i)*2*L)/(nx-1)
  rarray[i] = x
h5dset = h5grp.create_dataset('Data', data=rarray)

# 1d line in x
# Temperature boundary condition; NORTH
h5grp = h5file.create_group('x_line_node_centered_north')
h5grp.attrs['Dimension'] = np.string_('X')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [2.]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
nx = 2
rarray = np.zeros(nx,'=f8')
for i in range(nx):
  # T = 2 + x
  x = (float(i)*2*L)/(nx-1)
  rarray[i] = 2 + x
h5dset = h5grp.create_dataset('Data', data=rarray)

# 1d line in y
# Temperature boundary condition; EAST
h5grp = h5file.create_group('y_line_node_centered_east')
h5grp.attrs['Dimension'] = np.string_('Y')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [1.]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
ny = 2
rarray = np.zeros(ny,'=f8')
for j in range(ny):
  # T = 2 + 2*y
  y = (float(j)*L)/(ny-1)
  rarray[j] = 2 + 2*y
h5dset = h5grp.create_dataset('Data', data=rarray)