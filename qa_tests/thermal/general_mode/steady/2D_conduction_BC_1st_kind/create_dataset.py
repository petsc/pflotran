import sys
from h5py import *
import numpy as np

filename = 'dataset.h5'
h5file = File(filename,mode='w')

# 1d line in x
# Temperature boundary condition
h5grp = h5file.create_group('x_line_node_centered')
h5grp.attrs['Dimension'] = np.string_('X')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [1.]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
nx = 2
rarray = np.zeros(nx,'=f8')
for i in range(nx):
  rarray[i] = float(i)
h5dset = h5grp.create_dataset('Data', data=rarray)

# 1d line in y
# Temperature boundary condition
h5grp = h5file.create_group('y_line_node_centered')
h5grp.attrs['Dimension'] = np.string_('Y')
# Delta length between points [m]
h5grp.attrs['Discretization'] = [1.]
# Location of origin
h5grp.attrs['Origin'] = [0.]
# Load the dataset values
ny = 2
rarray = np.zeros(ny,'=f8')
for i in range(ny):
  rarray[i] = float(i)
h5dset = h5grp.create_dataset('Data', data=rarray)